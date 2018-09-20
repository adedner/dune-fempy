import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from math import pi,log10
from ufl import cos, sin, exp, as_vector, dx, grad, inner

# Crank Nicholson
theta = 0.5
deltaT = 0.01

from dune.fem import parameter
parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

order = 2
#grid = create.grid( "ALUSimplex", "../../data/hole2_larger.dgf", dimgrid=2 )

grid = create.grid("ALUCube",constructor=cartesianDomain([0,0],[3,1],[30,10]))
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

exact_u     = as_vector( [x[1] * (1.-x[1]), 0] )
exact_p     = as_vector( [ (-2*x[0] + 2)*0.1 ] )
# exact_u     = as_vector( [-1*cos(2.*pi*x[0])*sin(2.*pi*x[1])*exp(-8.*pi*pi*mu*time), sin(2.*pi*x[0])*cos(2.*pi*x[1])*exp(-8.*pi*pi*mu*time)])
# exact_p     = as_vector( [ -0.25*cos(4.*pi*x[0])*cos(4.*pi*x[1])*exp(-16.*pi*pi*mu*time) ] )
f           = as_vector( [0,]*grid.dimension )
f          += nu*exact_u
mainModel   = (deltaT*nu*dot(u,v) + deltaT*mu*inner(grad(u)+grad(u).T, grad(v)) - deltaT*dot(f,v)) * dx
gradModel   = -1*deltaT*inner( p[0]*Identity(grid.dimension), grad(v) ) * dx
divModel    = -div(u)*q[0] * dx
massModel   = inner(p,q) * dx
preconModel = inner(grad(p),grad(q)) * dx


# can use 'h1' or 'galerkin'
mainOp      = create.scheme("h1",spcU,(mainModel==0,DirichletBC(spcU,exact_u,1)), solver="cg")
# mainOp      = create.scheme("galerkin",(mainModel==0,DirichletBC(spcU,exact_u,1)))
gradOp      = create.operator("h1",gradModel)
divOp       = create.operator("galerkin",divModel)
massOp      = create.scheme("galerkin",massModel==0)
preconOp    = create.scheme("h1",preconModel==0)


mainOp.model.mu = 0.1
mainOp.model.nu = 1.0
# velocity = spcU.interpolate([0,]*spcU.dimRange, name="velocity")
# pressure = spcP.interpolate([0], name="pressure")
velocity = spcU.interpolate( lambda x: [ x[1] * ( 1.0 - x[ 1 ] ),0], name = "velocity", storage = "Istl" )
pressure = spcP.interpolate( lambda x: [(-2*x[0] + 2)*0.1], name = "pressure", storage = "Istl" )

grid.writeVTK("NavierStokes",
        pointdata={"pressure":pressure,
                   "exact_p":exact_p},
        pointvector={"velocity":velocity,
                     "exact_velo":exact_u},
        number=0
)

def compute():
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


    old_solution = velocity.copy();
    old_solution.name = "uOld"
    u_n = old_solution
    tau = 1
    solution = velocity, pressure

    # def plot(count):
    #     grid.writeVTK("NavierStokes",
    #             pointdata={"pressure":pressure, "rhsPress":rhsPress,
    #                        "exact_p":exact_p},
    #             pointvector={"velocity":velocity, "rhsVelo":rhsVelo,
    #                          "exact_velo":exact_u},
    #             number=count
    #     )


    def solveStokes(rhsVelo,r,d):
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
            velocity.axpy(-rho,xi)                #   velocity_.axpy(-rho,xi_);
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

    print( 'Solve step 1 - Stokes' )
    solveStokes(rhsVelo,r,d)
    print( 'Solve step 2 - Burgers' )
    solveBurger()
    print( 'Solve step 3 - Stokes' )
    solveStokes(rhsVelo,r,d)

    # plot()
    # vtk = grid.sequencedVTK("heat", pointdata=[velocity])
    # vtk()


endTime = 10
timeStep = deltaT
time = timeStep
counter = 1
while time < endTime:
    print( "Time is:", time )
    compute()
    grid.writeVTK("NavierStokes",
            pointdata={"pressure":pressure,
                       "exact_p":exact_p},
            pointvector={"velocity":velocity,
                         "exact_velo":exact_u},
            number=counter
    )
    # vtk = grid.writeVTK("ns_", pointdata=[velocity, pressure], number=counter )
    time += timeStep
    counter += 1
