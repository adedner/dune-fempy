import numpy
from scipy.sparse import bmat, linalg
import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, dot, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from dune.fem.operator import linear as linearOperator

order = 2
grid = create.grid("ALUCube",constructor=cartesianDomain([0,0],[3,1],[30,10]))
spcU = create.space("lagrange", grid, dimrange=grid.dimension, order=order, storage="fem")
spcP = create.space("lagrange", grid, dimrange=1, order=order-1, storage="fem")

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

# can also use 'operator' everywhere
mainOp      = create.scheme("galerkin",(mainModel==0,DirichletBC(spcU,exact_u,1)),spcU)
# mainOp      = create.scheme("h1",(mainModel==0,DirichletBC(spcU,exact_u,1)),spcU)
gradOp      = create.operator("h1",gradModel,spcP,spcU)
divOp       = create.operator("galerkin",divModel,spcU,spcP)
massOp      = create.scheme("galerkin",massModel==0,spcP)
preconOp    = create.scheme("h1",preconModel==0,spcP)

mainOp.model.mu = 0.1
mainOp.model.nu = 0.01

velocity = spcU.interpolate([0,]*spcU.dimRange, name="velocity")
pressure = spcP.interpolate([0], name="pressure")
rhsVelo  = velocity.copy()
rhsPress = pressure.copy()

sol_u  = velocity.as_numpy
sol_p  = pressure.as_numpy
rhs_u  = rhsVelo.as_numpy
rhs_p  = rhsPress.as_numpy
r      = numpy.copy(rhs_p)
d      = numpy.copy(rhs_p)
precon = numpy.copy(rhs_p)
xi     = numpy.copy(rhs_u)
A = linearOperator(mainOp).as_numpy
G = linearOperator(gradOp).as_numpy
D = linearOperator(divOp).as_numpy
M = linearOperator(massOp).as_numpy
P = linearOperator(preconOp).as_numpy
def Ainv(rhs,target): target[:] = linalg.spsolve(A,rhs)
def Minv(rhs,target): target[:] = linalg.spsolve(M,rhs)
def Pinv(rhs,target): target[:] = linalg.spsolve(P,rhs)

mainOp(velocity,rhsVelo)                # assembleRHS ( mainModel_, mainModel_.rightHandSide(), mainModel_.neumanBoundary(), rhsU_ );
rhs_u *= -1
xi[:] = G*sol_p                         # gradOperator_(pressure_, xi_);
rhs_u -= xi                             # rhs_ -= xi_;
mainOp.setConstraints(rhsVelo)          # mainOperator_.prepare( velocity_, rhsU_ );

def plot(count):
    grid.writeVTK("stokes",
            pointdata={"pressure":pressure, "rhsPress":rhsPress,
                       "exact_p":exact_p},
            pointvector={"velocity":velocity, "rhsVelo":rhsVelo,
                         "exact_velo":exact_u},
            number=count
    )


Ainv(rhs_u[:], sol_u[:])                # invMainOp( rhsU_, velocity_ );
rhs_p[:] = D*sol_u                      # divLinearOperator_( velocity_, rhsP_ );
Minv(rhs_p, r)                          # invMassOp( rhsP_, r_ );
if mainOp.model.nu > 0:                 # if (usePrecond_ && nu_>0.)
    precon.fill(0)                      #     precond_.clear();
    Pinv(rhs_p, precon)                 #   invPrecondOp( rhsP_, precond_ );
    r *= mainOp.model.mu                #   r_ *= mu_;
    r += mainOp.model.nu*precon         #   r_.axpy(nu_,precond_);
d[:] = r[:]                             # d_.assign(r_);
delta = numpy.dot(r,rhs_p)              # double delta = r_.scalarProductDofs(rhsP_);
assert delta >= 0                       # assert( delta >= 0 );
for m in range(100):                    # for (int m=0;m<100;++m)
    xi.fill(0)                          #   xi_.clear();
    rhs_u[:] = G*d                      #   gradLinearOperator_(d_, rhsU_);
    mainOp.setConstraints(\
       [0,]*grid.dimension, rhsVelo)    #   mainOperator_.prepare( mainModel_.zeroVelocity(), rhsU_ );
    Ainv(rhs_u[:], xi[:])               #   invMainOp( rhsU_, xi_ );
    rhs_p[:] = D*xi                     #   divLinearOperator_( xi_, rhsP_ );
    rho = delta / numpy.dot(d,rhs_p)    #   double rho = delta / d_.scalarProductDofs(rhsP_);
    sol_p += rho*d                      #   pressure_.axpy(rho,d_);
    sol_u -= rho*xi                     #   velocity_.axpy(-rho,xi_);
    rhs_p[:] = D*sol_u                  #   divLinearOperator_( velocity_, rhsP_ );
    Minv(rhs_p[:],r[:])                 #   invMassOp( rhsP_, r_ );
    if mainOp.model.nu > 0:             #   if (usePrecond_ && nu_>0.)
        precon.fill(0)                  #     precond_.clear();
        Pinv(rhs_p,precon)              #     invPrecondOp( rhsP_, precond_ );
        r *= mainOp.model.mu            #     r_ *= mu_;
        r += mainOp.model.nu*precon     #     r_.axpy(nu_,precond_);
    oldDelta = delta                    #     double oldDelta = delta;
    delta = numpy.dot(r,rhs_p)          #     delta = r_.scalarProductDofs(rhsP_);
    print("delta:",delta)               #     std::cout << "delta: " << delta << std::endl;
    if delta < 1e-14: break             #     if ( delta < solverEps_*10. ) break;
    gamma = delta/oldDelta              #     double gamma = delta/oldDelta;
    d *= gamma                          #     d_ *= gamma;
    d += r                              #     d_ += r_;
plot(0)
