# NOTE: there is some issue with failing convergence when using solver=cg -
# it should work...
import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, dot, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from dune.fem import parameter
from dune.fem.space import lagrange  as lagrangeSpace
from dune.fem.scheme import h1 as h1Scheme
from dune.fem.operator import galerkin as galerkinOperator
from dune.fem.operator import h1 as h1Operator
from dune.fem.operator import linear as linearOperator

parameter.append({"fem.verboserank": -1})

storage="petsc"
# storage="istl"
# storage="fem"
order = 2
grid = create.grid("ALUCube",constructor=cartesianDomain([0,0],[3,1],[30,10]))

spcU = lagrangeSpace(grid, dimrange=grid.dimension, order=order, storage=storage)
spcP = lagrangeSpace(grid, dimrange=1, order=order-1, storage=storage)

cell  = spcU.cell()
x     = SpatialCoordinate(cell)
mu    = NamedConstant(cell, "mu")
nu    = NamedConstant(cell, "nu")
u     = TrialFunction(spcU)
v     = TestFunction(spcU)
p     = TrialFunction(spcP)
q     = TestFunction(spcP)
exact_u     = as_vector( [x[1] * (1.-x[1]), 0] )
exact_p     = as_vector( [ (-2*x[0] + 2)*mu ] )
f           = as_vector( [0,]*grid.dimension )
f          += nu*exact_u
mainModel   = (nu*dot(u,v) + mu*inner(grad(u)+grad(u).T, grad(v)) - dot(f,v)) * dx
gradModel   = -inner( p[0]*Identity(grid.dimension), grad(v) ) * dx
divModel    = -div(u)*q[0] * dx
massModel   = inner(p,q) * dx
preconModel = inner(grad(p),grad(q)) * dx

solverParameters = {}
if storage == "istl" :
    solverParameters = {"istl.preconditioning.method": "ilu", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2}

# can use 'h1' or 'galerkin'
mainOp      = create.scheme("galerkin",(mainModel==0,DirichletBC(spcU,exact_u,1)),
# mainOp      = h1Scheme((mainModel==0,DirichletBC(spcU,exact_u,1),spcU),
                        solver="cg", parameters=solverParameters)
gradOp      = h1Operator(gradModel)
divOp       = galerkinOperator(divModel)
massOp      = h1Scheme(massModel==0,
                       solver="gmres", parameters=solverParameters)
preconOp    = h1Scheme(preconModel==0,
                       solver="cg", parameters=solverParameters)

mainOp.model.mu = 0.1
mainOp.model.nu = 0.01

velocity = spcU.interpolate([0,]*spcU.dimRange, name="velocity")
pressure = spcP.interpolate([0], name="pressure")
rhsVelo  = velocity.copy()
rhsPress = pressure.copy()

r      = rhsPress.copy()
d      = rhsPress.copy()
precon = rhsPress.copy()
xi     = rhsVelo.copy()

print("ASSEMBLE all matrices")
A = linearOperator(mainOp)
G = linearOperator(gradOp)
D = linearOperator(divOp)
M = linearOperator(massOp)
P = linearOperator(preconOp)

solver = {"method":"gmres","verbose":True,
          "relativetol":1e-10,"absolutetol":1e-10}
Ainv   = mainOp.inverseLinearOperator(A,parameters=solver)
Minv   = massOp.inverseLinearOperator(M,solver)
Pinv   = preconOp.inverseLinearOperator(P,solver)

def plot(count=None):
    grid.writeVTK("Stokes",
            pointdata={"pressure":pressure, "rhsPress":rhsPress,
                       "exact_p":exact_p},
            pointvector={"velocity":velocity, "rhsVelo":rhsVelo,
                         "exact_velo":exact_u},
            number=count
    )

                                          # dune-fem-howto version:
                                          #############################
mainOp(velocity,rhsVelo)                  # assembleRHS ( mainModel_, mainModel_.rightHandSide(), mainModel_.neumanBoundary(), rhsU_ );
rhsVelo *= -1
G(pressure,xi)                            # gradOperator_(pressure_, xi_);
rhsVelo -= xi                             # rhs_ -= xi_;
mainOp.setConstraints(rhsVelo)            # mainOperator_.prepare( velocity_, rhsU_ );

Ainv(rhsVelo, velocity)                   # invMainOp( rhsU_, velocity_ );
D(velocity,rhsPress)                      # divLinearOperator_( velocity_, rhsP_ );
Minv(rhsPress, r)                         # invMassOp( rhsP_, r_ );

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
    if delta < 1e-9: break                #     if ( delta < solverEps_*10. ) break;
    gamma = delta/oldDelta                #     double gamma = delta/oldDelta;
    d *= gamma                            #     d_ *= gamma;
    d += r                                #     d_ += r_;
plot()

# add a version that uses methods from petsc4py
# if storage == "petsc" ...
