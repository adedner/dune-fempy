import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from math import pi,log10
from dune.fem.function import integrate

from ufl import cos, sin, exp, as_vector, dx, grad, inner,sqrt

from dune.fem import parameter
parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

endTime = 0.1
order   = 2

grid = create.grid("ALUCube",constructor=cartesianDomain([-1,-1],[1,1],[50,50]))
spcU = create.space("lagrange", grid, dimrange=grid.dimension, order=order, storage="istl")
spcP = create.space("lagrange", grid, dimrange=1, order=order-1, storage="istl")
x     = SpatialCoordinate(spcU)
dt    = NamedConstant(spcU, "dt")
time  = NamedConstant(spcU, "t")     # current time
u     = TrialFunction(spcU)
v     = TestFunction(spcU)
p     = TrialFunction(spcP)
q     = TestFunction(spcP)

gamma_t = lambda t: exp(-2.*pi*pi*mu*t)
exact_u = as_vector( [-1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t(time), sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*mu*time)])
exact_p = as_vector( [ -0.25*(cos(2.*pi*x[0])+cos(2.*pi*x[1]))*exp(-4.*pi*pi*mu*time) ] )

exact_uend = as_vector( [-1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*exp(-2.*pi*pi*mu*0.1), sin(1.*pi*x[0])*cos(1.*pi*x[1])*gamma_t(endTime)

f           = as_vector( [2.*pi*pi*mu*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t(time) , -2.*pi*pi*mu*sin(1.*pi*x[0])*cos(1.*pi*x[1])*gamma_t(time) ] )
f          += as_vector( [0.25*2.*pi*sin(2.*pi*x[0])*exp(-4.*pi*pi*mu*time),  0.25*2.*pi*sin(2.*pi*x[1])*exp(-4.*pi*pi*mu*time)] )
f          -= as_vector( [mu*2*gamma_t(time)*pi*pi*cos(1.*pi*x[0])*sin(1.*pi*x[1]), mu*-2*gamma_t(time)*pi*pi*sin(1.*pi*x[0])*cos(1.*pi*x[1])] )
f          += as_vector( [cos(1.*pi*x[0])*sin(1.*pi*x[0])*(cos(1.*pi*x[1])*cos(1.*pi*x[1])-sin(1.*pi*x[1])*sin(1.*pi*x[1])),cos(1.*pi*x[1])*sin(1.*pi*x[1])*(cos(1.*pi*x[0])*cos(1.*pi*x[0])-sin(1.*pi*x[0])*sin(1.*pi*x[0]))] )
f          += nu*exact_u

##############################################

class Stokes:
    def __init__(self,spc,forcing):
        spcU = spc[0]
        spcP = spc[1]
        mu = NamedConstant(spcU, "mu")
        nu = NamedConstant(spcU, "nu")
        self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")
        mainModel    = (mu*dot(u,v) + nu*inner(grad(u)+grad(u).T, grad(v))) * dx
        gradModel    = -1*inner( p[0]*Identity(grid.dimension), grad(v) ) * dx
        divModel     = -1*inner(div(u),q[0]) * dx
        massModel    = inner(p,q) * dx
        preconModel  = inner(grad(p),grad(q)) * dx
        f = (dot(forcing,v) +
             mu*dot(oldVelo,v) - nu*inner(grad(oldVelo)+grad(oldVelo).T, grad(v)) -
             inner(grad(oldVelo), outer(v, oldVelo))) * dx

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
        mainOp(velocity,rhsVelo)                  # assembleRHS ( mainModel_, mainModel_.rightHandSide(), mainModel_.neumanBoundary(), rhsU_ );
        rhsVelo *= -1
        G(pressure,xi)                            # gradOperator_(pressure_, xi_);
        rhsVelo -= xi                             # rhs_ -= xi_;
        mainOp.setConstraints(rhsVelo)            # mainOperator_.prepare( velocity_, rhsU_ );

        Ainv(rhsVelo, velocity)                   # invMainOp( rhsU_, velocity_ );
        D(velocity,rhsPress)                      # divLinearOperator_( velocity_, rhsP_ );
        Minv(rhsPress, r)                         # invMassOp( rhsP_, r_ );

        # steps = int(0.1 / deltaT)
        # for n in range(1,10):
        if mainOp.model.nu > 0:                   # if (usePrecond_ && nu_>0.)
            precon.clear()                        #   precond_.clear();
            Pinv(rhsPress, precon)                #   invPrecondOp( rhsP_, precond_ );
            r *= mainOp.model.mu                  #   r_ *= mu_;
            r.axpy(mainOp.model.nu,precon)        #   r_.axpy(nu_,precond_);
        d.assign(r)                               # d_.assign(r_);
        delta = r.scalarProductDofs(rhsPress)     # double delta = r_.scalarProductDofs(rhsP_);
        print("delta:",delta,flush=True)          # std::cout << "delta: " << delta << std::endl;
        assert delta >= 0                         # assert( delta >= 0 );
        for m in range(100):                      # for (int m=0;m<100;++m)
            xi.clear()                            #   xi_.clear();
            G(d,rhsVelo)                          #   gradLinearOperator_(d_, rhsU_);
            mainOp.setConstraints(\
                [0,]*grid.dimension, rhsVelo)     #   mainOperator_.prepare( mainModel_.zeroVelocity(), rhsU_ );
            Ainv(rhsVelo, xi)                     #   invMainOp( rhsU_, xi_ );
            D(xi,rhsPress)                        #   divLinearOperator_( xi_, rhsP_ );
            rho = delta /\
               d.scalarProductDofs(rhsPress)      #   double rho = delta / d_.scalarProductDofs(rhsP_);
            pressure.axpy(rho,d)                  #   pressure_.axpy(rho,d_);
            velocity.axpy(-1*rho,xi)                #   velocity_.axpy(-rho,xi_);
            D(velocity, rhsPress)                 #   divLinearOperator_( velocity_, rhsP_ );
            Minv(rhsPress,r)                      #   invMassOp( rhsP_, r_ );
            if mainOp.model.nu > 0:               #   if (usePrecond_ && nu_>0.)
                precon.clear()                    #     precond_.clear();
                Pinv(rhsPress,precon)             #     invPrecondOp( rhsP_, precond_ );
                r *= mainOp.model.mu              #     r_ *= mu_;
                r.axpy(mainOp.model.nu,precon)    #     r_.axpy(nu_,precond_);
            oldDelta = delta                      #     double oldDelta = delta;
            delta = r.scalarProductDofs(rhsPress) #     delta = r_.scalarProductDofs(rhsP_);
            print("delta:",delta,flush=True)      #     std::cout << "delta: " << delta << std::endl;
            if delta < 1e-14: break               #     if ( delta < solverEps_*10. ) break;
            gamma = delta/oldDelta                #     double gamma = delta/oldDelta;
            d *= gamma                            #     d_ *= gamma;
            d += r                                #     d_ += r_;


    # a = (inner(u - u_n, v) + inner(grad(u), outer(v, u_n))) * dx
    a = (inner(u - u_n, v) + deltaT * inner(grad(u), outer(v, u))) * dx
    model = create.model("integrands", grid, a == 0)

    solverParameter={"fem.solver.newton.linabstol": 1e-11,
                     "fem.solver.newton.linreduction": 1e-11,
                     "fem.solver.newton.tolerance": 1e-10,
                     "fem.solver.newton.verbose": "true",
                     "fem.solver.newton.linear.verbose": "false"}

    scheme = create.scheme("galerkin", a == 0, spcU, solver="gmres", parameters = solverParameter)

    def solveBurger():
        old_solution.assign(velocity)
        scheme.solve(target=velocity)
        # vtk()

    endTime = 0.1
    timeStep = deltaT
    time = timeStep
    counter = 0
    while time < endTime:
        print( "Time is:", time )
        print( 'Solve step 1 - Stokes' )
        solveStokes()
        # print( 'Solve step 2 - Burgers' )
        solveBurger()
        # print( 'Solve step 3 - Stokes' )
        solveStokes()

        l2error_fn = dot(velocity - exact_uend, velocity - exact_uend)
        l2error = sqrt( integrate(grid, l2error_fn, 5)[0] )
        print('|u_h - u| =', l2error)


        vtk()


        time += timeStep


velocity = spcU.interpolate(exact_u, name="velocity")
pressure = spcP.interpolate(exact_p, name="pressure")
# velocity = spcU.interpolate( exact_u, name = "velocity", storage = "Istl" )
# pressure = spcP.interpolate( exact_p, name = "pressure", storage = "Istl" )
vtk = grid.sequencedVTK("navierfracstep_model2",pointdata={"pressure":pressure,
           "exact_p":exact_p}, pointvector={"velocity":velocity})
vtk()

# old_solution = velocity.copy();
# old_solution.name = "uOld"
# u_n = old_solution

        old_solution = velocity.copy();
        old_solution.name = "uOld"
        u_n = old_solution
        # tau = 1
        # solution = velocity, pressure

compute()
