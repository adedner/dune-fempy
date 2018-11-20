from dune import create
from dune.fem.function import integrate


class Stokes:
    def __init__(self,spc,forcing):
        spcU = spc[0]
        spcP = spc[1]
        u  = TrialFunction(spcU)
        v  = TestFunction(spcU)
        p  = TrialFunction(spcP)
        q  = TestFunction(spcP)
        mu = NamedConstant(spcU, "mu")
        nu = NamedConstant(spcU, "nu")
        self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")
        mainModel    = (mu*dot(u,v) + nu*inner(grad(u)+grad(u).T, grad(v))) * dx
        gradModel    = -1*inner( p[0]*Identity(grid.dimension), grad(v) ) * dx
        divModel     = -1*inner(div(u),q[0]) * dx
        massModel    = inner(p,q) * dx
        preconModel  = inner(grad(p),grad(q)) * dx
        f = (dot(forcing,v) +
             mu*dot(self.oldVelo,v) - nu*inner(grad(self.oldVelo)+grad(self.oldVelo).T, grad(v)) -
             inner(grad(self.oldVelo), outer(v, self.oldVelo))) * dx

        self.mainOp   = create.scheme("galerkin",spcU,(mainModel==f,DirichletBC(spcU,exact_u,1)), solver="cg")
        self.gradOp   = create.operator("galerkin",gradModel)
        self.divOp    = create.operator("galerkin",divModel)
        self.massOp   = create.scheme("galerkin",massModel==0)
        self.preconOp = create.scheme("galerkin",preconModel==0)

        self.rhsVelo  = spcU.interpolate([0,0],name="rhsU") # TODO: extract correct dimension
        self.rhsPress = spcP.interpolate([0],name="rhsP")
        self.r        = rhsPress.copy()
        self.d        = rhsPress.copy()
        self.precon   = rhsPress.copy()
        self.xi       = rhsVelo.copy()

    def prepare(self,t,mu,nu,velo):
        try:
            self.mainOp.model.time = t
        except: pass
        self.mainOp.model.mu = mu
        self.mainOp.model.nu = nu
        self.oldVelo.assign(velo)
        self.A = self.mainOp.assemble(self.rhsVelo)
        self.G = self.gradOp.assemble(self.rhsPress)
        self.D = self.divOp.assemble(self.rhsVelo)
        self.M = self.massOp.assemble(self.rhsPress)
        self.P = self.preconOp.assemble(self.rhsPress)
        solver = {"krylovmethod":"cg","fem.solver.verbose":0}
        self.Ainv   = self.mainOp.inverseLinearOperator(A,1e-10,parameters=solver)
        self.Minv   = self.massOp.inverseLinearOperator(M,1e-10,solver)
        self.Pinv   = self.preconOp.inverseLinearOperator(P,1e-10,solver)

    def solve(self,target):
        velocity = target[0]
        pressure = target[1]
        self.mainOp(velocity,self.rhsVelo)                  # assembleRHS ( mainModel_, mainModel_.rightHandSide(), mainModel_.neumanBoundary(), rhsU_ );
        self.rhsVelo *= -1
        self.G(pressure,xi)                            # gradOperator_(pressure_, xi_);
        self.rhsVelo -= xi                             # rhs_ -= xi_;
        self.mainOp.setConstraints(self.rhsVelo)            # mainOperator_.prepare( velocity_, rhsU_ );

        self.Ainv(self.rhsVelo, self.velocity)                   # invMainOp( rhsU_, velocity_ );
        self.D(velocity,self.rhsPress)                      # divLinearOperator_( velocity_, rhsP_ );
        self.Minv(self.rhsPress, self.r)                         # invMassOp( rhsP_, r_ );

        # steps = int(0.1 / deltaT)
        # for n in range(1,10):
        if self.mainOp.model.mu > 0:                   # if (usePrecond_ && nu_>0.)
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
            self.mainOp.setConstraints(\
                [0,]*grid.dimension, self.rhsVelo)     #   mainOperator_.prepare( mainModel_.zeroVelocity(), rhsU_ );
            self.Ainv(self.rhsVelo, self.xi)                     #   invMainOp( rhsU_, xi_ );
            self.D(self.xi,self.rhsPress)                        #   divLinearOperator_( xi_, rhsP_ );
            self.rho = delta /\
               self.d.scalarProductDofs(self.rhsPress)      #   double rho = delta / d_.scalarProductDofs(rhsP_);
            pressure.axpy(self.rho,self.d)                  #   pressure_.axpy(rho,d_);
            velocity.axpy(-1*self.rho,self.xi)                #   velocity_.axpy(-rho,xi_);
            self.D(velocity, self.rhsPress)                 #   divLinearOperator_( velocity_, rhsP_ );
            self.Minv(self.rhsPress,self.r)                      #   invMassOp( rhsP_, r_ );
            if self.mainOp.model.mu > 0:               #   if (usePrecond_ && nu_>0.)
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
            self.d += r                                #     d_ += r_;
