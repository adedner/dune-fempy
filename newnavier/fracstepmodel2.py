import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from math import pi,log10
from dune.fem.function import integrate

from ufl import cos, sin, exp, as_vector, dx, grad, inner,sqrt

deltaT = 0.001

from dune.fem import parameter
parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "jacobi", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

order = 2
# grid = create.grid( "ALUSimplex", "../../data/unitcube-2d.dgf", dimgrid=2 )
# grid.hierarchicalGrid.globalRefine(2)

grid = create.grid("ALUCube",constructor=cartesianDomain([-1,-1],[1,1],[100,100]))
# grid = create.grid("ALUCube",constructor=cartesianDomain([0,0],[1,1],[10,10]))
spcU = create.space("lagrange", grid, dimrange=grid.dimension, order=order, storage="istl")
spcP = create.space("lagrange", grid, dimrange=1, order=order-1, storage="istl")
cell  = spcU.cell()
x     = SpatialCoordinate(cell)
mu    = NamedConstant(cell, "mu")
nu    = NamedConstant(cell, "nu")
u     = TrialFunction(spcU)
v     = TestFunction(spcU)
p     = TrialFunction(spcP)
q     = TestFunction(spcP)
time     = NamedConstant(triangle, "t")     # current time


alpha  = 0.5
beta   = 0.5
theta1 = 1/3
theta2 = 1/3
mu     = 0.1
nu     = 0.0

# exact_u     = as_vector( [x[1] * (1.-x[1]), 0] )
# exact_p     = as_vector( [ (-2*x[0] + 2)*0.1 ] )
gamma_t = exp(-2.*pi*pi*mu*time)
exact_u     = as_vector( [-1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t, sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*mu*time)])
exact_p     = as_vector( [ -0.25*(cos(2.*pi*x[0])+cos(2.*pi*x[1]))*exp(-4.*pi*pi*mu*time) ] )

exact_uend     = as_vector( [-1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*exp(-2.*pi*pi*mu*0.1), sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*mu*0.1)])


f           = as_vector( [2.*pi*pi*mu*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t , -2.*pi*pi*mu*sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*mu*time) ] )
f           += as_vector( [0.25*2.*pi*sin(2.*pi*x[0])*exp(-4.*pi*pi*mu*time),  0.25*2.*pi*sin(2.*pi*x[1])*exp(-4.*pi*pi*mu*time)] )

f           -= as_vector( [mu*2*gamma_t*pi*pi*cos(1.*pi*x[0])*sin(1.*pi*x[1]), mu*-2*gamma_t*pi*pi*sin(1.*pi*x[0])*cos(1.*pi*x[1])] )
f           += as_vector( [cos(1.*pi*x[0])*sin(1.*pi*x[0])*(cos(1.*pi*x[1])*cos(1.*pi*x[1])-sin(1.*pi*x[1])*sin(1.*pi*x[1])),cos(1.*pi*x[1])*sin(1.*pi*x[1])*(cos(1.*pi*x[0])*cos(1.*pi*x[0])-sin(1.*pi*x[0])*sin(1.*pi*x[0]))] )

f          += nu*exact_u


velocity = spcU.interpolate(exact_u, name="velocity")
pressure = spcP.interpolate(exact_p, name="pressure")
# velocity = spcU.interpolate( exact_u, name = "velocity", storage = "Istl" )
# pressure = spcP.interpolate( exact_p, name = "pressure", storage = "Istl" )
vtk = grid.sequencedVTK("visc1navierfracstep_model2",pointdata={"pressure":pressure,
           "exact_p":exact_p}, pointvector={"velocity":velocity})
vtk()

# old_solution = velocity.copy();
# old_solution.name = "uOld"
# u_n = old_solution

def compute():
    # mainOp.model.alpha = 0.5
    # mainOp.model.beta = 0.5
    # mainOp.model.theta1 = 1/3
    # mainOp.model.theta2 = 1/3

    def solveStokes():
        old_solution = velocity.copy();
        old_pressure = pressure.copy();

        old_solution.name = "uOld"
        u_n = old_solution
        old_pressure.name = "pOld"
        p_n = old_pressure
        # tau = 1
        # solution = velocity, pressure
        mainModel   =  (inner(u-u_n,v) + theta1 * deltaT * alpha * mu * inner(grad(u)+grad(u).T, grad(v))) * dx
        mainModel   += (theta1 * deltaT * beta * mu * inner(grad(u_n)+grad(u_n).T, grad(v)) - theta1 * deltaT * dot(f,v)) * dx
        mainModel   += theta1 * deltaT * inner(grad(u_n), outer(v, u_n)) * dx
        gradModel   =  -1* theta1 *deltaT * inner( p[0] * Identity(grid.dimension), grad(v) ) * dx
        divModel    =  -1 * inner(div(u),q[0]) * dx
        massModel   =  inner(p,q) * dx
        preconModel =  inner(grad(p),grad(q)) * dx


        # can use 'h1' or 'galerkin'
        mainOp      = create.scheme("h1",spcU,(mainModel==0,DirichletBC(spcU,exact_u,1)), solver="cg")
        # mainOp      = create.scheme("galerkin",(mainModel==0,DirichletBC(spcU,exact_u,1)))
        gradOp      = create.operator("h1",gradModel)
        divOp       = create.operator("galerkin",divModel)
        massOp      = create.scheme("galerkin",massModel==0)
        preconOp    = create.scheme("h1",preconModel==0)


        # mainOp.model.mu = 0.01
        # mainOp.model.nu = 0.0

        rhsVelo  = velocity.copy()
        rhsPress = pressure.copy()

        r      = rhsPress.copy()
        d      = rhsPress.copy()
        precon = rhsPress.copy()
        xi     = rhsVelo.copy()

        # Question: should assemble method also provide the affine shift?
        A      = mainOp.assemble(velocity)
        G      = gradOp.assemble(pressure)
        D      = divOp.assemble(velocity)
        M      = massOp.assemble(pressure)
        P      = preconOp.assemble(pressure)
        solver = {"krylovmethod":"cg","fem.solver.verbose":0}
        Ainv   = mainOp.inverseLinearOperator(A,1e-10,parameters=solver)
        Minv   = massOp.inverseLinearOperator(M,1e-10,solver)
        Pinv   = preconOp.inverseLinearOperator(P,1e-10,solver)


        # old_solution = velocity.copy();
        # old_solution.name = "uOld"
        # u_n = old_solution
        tau = 1
        solution = velocity, pressure

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
        if nu > 0:                   # if (usePrecond_ && nu_>0.)
            precon.clear()                        #   precond_.clear();
            Pinv(rhsPress, precon)                #   invPrecondOp( rhsP_, precond_ );
            r *= mu                  #   r_ *= mu_;
            r.axpy(nu,precon)        #   r_.axpy(nu_,precond_);
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
            if nu > 0:               #   if (usePrecond_ && nu_>0.)
                precon.clear()                    #     precond_.clear();
                Pinv(rhsPress,precon)             #     invPrecondOp( rhsP_, precond_ );
                r *= mu              #     r_ *= mu_;
                r.axpy(nu,precon)    #     r_.axpy(nu_,precond_);
            oldDelta = delta                      #     double oldDelta = delta;
            delta = r.scalarProductDofs(rhsPress) #     delta = r_.scalarProductDofs(rhsP_);
            print("delta:",delta,flush=True)      #     std::cout << "delta: " << delta << std::endl;
            if delta < 1e-14: break               #     if ( delta < solverEps_*10. ) break;
            gamma = delta/oldDelta                #     double gamma = delta/oldDelta;
            d *= gamma                            #     d_ *= gamma;
            d += r                                #     d_ += r_;



    def solveBurger():
        old_solution = velocity.copy();
        old_pressure = pressure.copy();

        old_solution.name = "uOld"
        u_n = old_solution
        old_pressure.name = "pOld"
        p_n = old_pressure
        # a = (inner(u - u_n, v) + inner(grad(u), outer(v, u_n))) * dx
        a = (inner(u - u_n, v) + theta2 * deltaT * beta * mu * inner(grad(u)+grad(u).T, grad(v)) + theta2 * deltaT * inner(grad(u), outer(v, u))) * dx
        a += ( theta2 * deltaT * alpha * mu * inner(grad(u_n)+grad(u_n).T, grad(v))) * dx
        a += (-1* theta2 *deltaT * inner(p_n[0],div(v)) )* dx


        model = create.model("integrands", grid, a == 0)

        solverParameter={"fem.solver.newton.linabstol": 1e-11,
                         "fem.solver.newton.linreduction": 1e-11,
                         "fem.solver.newton.tolerance": 1e-10,
                         "fem.solver.newton.verbose": "true",
                         "fem.solver.newton.linear.verbose": "false"}

        scheme = create.scheme("galerkin", a == 0, spcU, solver="gmres", parameters = solverParameter)

        old_solution.assign(velocity)
        old_pressure.assign(pressure)

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
        time += (1/3) * timeStep
        print( 'Solve step 2 - Burgers' )
        solveBurger()
        time += (1/3) * timeStep
        print( 'Solve step 3 - Stokes' )
        solveStokes()

        l2error_fn = dot(velocity - exact_uend, velocity - exact_uend)
        l2error = sqrt( integrate(grid, l2error_fn, 5)[0] )
        print('|u_h - u| =', l2error)
        vtk()
        time += (1/3) * timeStep

compute()
