from dune import create
from dune.fem.function import integrate
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant

class Stokes:
    def __init__(self,spc,forcing,bcs,Re,withBurgers=False):
        spcU = spc[0]
        spcP = spc[1]
        # TODO add dimDomain property to space
        dimension = spcU.grid.dimension
        u  = TrialFunction(spcU)
        v  = TestFunction(spcU)
        p  = TrialFunction(spcP)
        q  = TestFunction(spcP)
        mu = NamedConstant(spcU, "mu")
        nu = NamedConstant(spcU, "nu")
        x     = SpatialCoordinate(spcU)
        time  = NamedConstant(spcU, "time")     # current time
        self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")

        mainModel    = (mu*dot(u,v) + nu/Re*inner(grad(u)+grad(u).T, grad(v))) * dx
        gradModel    = -1*inner( p[0]*Identity(dimension), grad(v) ) * dx
        divModel     = -1*inner(div(u),q[0]) * dx
        massModel    = inner(p,q) * dx
        preconModel  = inner(grad(p),grad(q)) * dx
        if forcing is not None:
            f = dot(forcing,v) * dx
        else:
            f = 0
        if withBurgers:
             # self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")
             f += (mu*dot(self.oldVelo,v) - (1-nu)/Re*inner(grad(self.oldVelo)+grad(self.oldVelo).T, grad(v)) -
             inner(grad(self.oldVelo), outer(v, self.oldVelo))) * dx

        # TODO give better diagonostics when no equation is passed to scheme
        # TODO galerkin doesn't work here
        self.mainOp   = create.scheme("galerkin",[mainModel==f,*bcs],spcU)
        # self.mainOp   = create.scheme("h1",spcU,[mainModel==f,*bcs], solver="cg")
        self.gradOp   = create.operator("galerkin",gradModel)
        self.divOp    = create.operator("galerkin",divModel)
        self.massOp   = create.scheme("galerkin",massModel==0)
        self.preconOp = create.scheme("galerkin",preconModel==0)

        self.rhsVelo  = spcU.interpolate([0,]*spcU.dimRange,name="rhsU")
        self.rhsPress = spcP.interpolate([0],name="rhsP")
        self.r        = self.rhsPress.copy()
        self.d        = self.rhsPress.copy()
        self.precon   = self.rhsPress.copy()
        self.xi       = self.rhsVelo.copy()

        self.G = self.gradOp.assemble(self.rhsPress)
        self.D = self.divOp.assemble(self.rhsVelo)
        self.M = self.massOp.assemble(self.rhsPress)
        self.P = self.preconOp.assemble(self.rhsPress)
        self.solver = {"krylovmethod":"cg","fem.solver.verbose":0}
        self.Minv   = self.massOp.inverseLinearOperator(self.M,1e-10,self.solver)
        self.Pinv   = self.preconOp.inverseLinearOperator(self.P,1e-10,self.solver)

    def prepare(self,time=0,mu=0,nu=1,oldSolution=None):
        try:
            self.mainOp.model.time = time
        except: pass
        self.mainOp.model.mu = mu
        self.mainOp.model.nu = nu
        if oldSolution:
            self.oldVelo=oldSolution[0].copy()
            # self.oldVelo.assign(oldSolution[0])
        self.A = self.mainOp.assemble(self.rhsVelo)
        self.Ainv   = self.mainOp.inverseLinearOperator(self.A,1e-10,parameters=self.solver)

    def solve(self,target):
        velocity = target[0]
        pressure = target[1]
        self.mainOp(velocity,self.rhsVelo)                  # assembleRHS ( mainModel_, mainModel_.rightHandSide(), mainModel_.neumanBoundary(), rhsU_ );
        self.rhsVelo *= -1
        self.G(pressure,self.xi)                            # gradOperator_(pressure_, xi_);
        self.rhsVelo -= self.xi                             # rhs_ -= xi_;
        # try:
        self.mainOp.setConstraints(self.rhsVelo)        # mainOperator_.prepare( velocity_, rhsU_ );
        # except: pass
        self.Ainv(self.rhsVelo, velocity)                   # invMainOp( rhsU_, velocity_ );
        self.D(velocity,self.rhsPress)                      # divLinearOperator_( velocity_, rhsP_ );
        self.Minv(self.rhsPress, self.r)                         # invMassOp( rhsP_, r_ );

        # steps = int(0.1 / deltaT)
        # for n in range(1,10):
        if False: # self.mainOp.model.mu > 0:                   # if (usePrecond_ && nu_>0.)
            self.precon.clear()                        #   precond_.clear();
            self.Pinv(self.rhsPress, self.precon)                #   invPrecondOp( rhsP_, precond_ );
            self.r *= self.mainOp.model.nu                  #   r_ *= mu_;
            self.r.axpy(self.mainOp.model.mu,self.precon)        #   r_.axpy(nu_,precond_);
        self.d.assign(self.r)                               # d_.assign(r_);
        delta = self.r.scalarProductDofs(self.rhsPress)     # double delta = r_.scalarProductDofs(rhsP_);
        print("delta:",delta,flush=True)          # std::cout << "delta: " << delta << std::endl;
        assert delta >= 0                         # assert( delta >= 0 );
        for m in range(100):                      # for (int m=0;m<100;++m)
            self.xi.clear()                            #   xi_.clear();
            self.G(self.d,self.rhsVelo)                          #   gradLinearOperator_(d_, rhsU_);
            # try:
            self.mainOp.setConstraints(\
                velocity, self.rhsVelo)     #   mainOperator_.prepare( mainModel_.zeroVelocity(), rhsU_ );
            # except: pass
            self.Ainv(self.rhsVelo, self.xi)                     #   invMainOp( rhsU_, xi_ );
            self.D(self.xi,self.rhsPress)                        #   divLinearOperator_( xi_, rhsP_ );
            self.rho = delta /\
               self.d.scalarProductDofs(self.rhsPress)      #   double rho = delta / d_.scalarProductDofs(rhsP_);
            pressure.axpy(self.rho,self.d)                  #   pressure_.axpy(rho,d_);
            velocity.axpy(-1*self.rho,self.xi)                #   velocity_.axpy(-rho,xi_);
            self.D(velocity, self.rhsPress)                 #   divLinearOperator_( velocity_, rhsP_ );
            self.Minv(self.rhsPress,self.r)                      #   invMassOp( rhsP_, r_ );
            if False: # self.mainOp.model.mu > 0:               #   if (usePrecond_ && nu_>0.)
                self.precon.clear()                    #     precond_.clear();
                self.Pinv(self.rhsPress,self.precon)             #     invPrecondOp( rhsP_, precond_ );
                self.r *= self.mainOp.model.nu              #     r_ *= mu_;
                self.r.axpy(self.mainOp.model.mu,self.precon)    #     r_.axpy(nu_,precond_);
            oldDelta = delta                      #     double oldDelta = delta;
            delta = self.r.scalarProductDofs(self.rhsPress) #     delta = r_.scalarProductDofs(rhsP_);
            print("delta:",delta,flush=True)      #     std::cout << "delta: " << delta << std::endl;
            if delta < 1e-14: break               #     if ( delta < solverEps_*10. ) break;
            gamma = delta/oldDelta                #     double gamma = delta/oldDelta;
            self.d *= gamma                            #     d_ *= gamma;
            self.d += self.r                                #     d_ += r_;
