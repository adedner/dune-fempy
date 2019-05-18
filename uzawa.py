# NOTE: there is some issue with failing convergence when using solver=cg -
# it should work...
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
from dune.grid import structuredGrid, cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, dot, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import Constant, DirichletBC
from dune.fem.space import lagrange  as lagrangeSpace
from dune.fem.operator import galerkin as galerkinOperator
from dune.fem.operator import linear as linearOperator
from dune.fem.scheme import galerkin as galerkinScheme

from dune.fem import parameter
parameter.append({"fem.verboserank": 0})


order = 2
# Note: structuredGrid fails in precon step!
# grid = structuredGrid([0,0],[3,1],[30,10])
from dune.alugrid import aluCubeGrid as leafGridView
grid = leafGridView( cartesianDomain([0,0],[3,1],[30,10]) )

spcU = lagrangeSpace(grid, dimRange=grid.dimension, order=order, storage="petsc")
spcP = lagrangeSpace(grid, dimRange=1, order=order-1, storage="petsc")

cell  = spcU.cell()
x     = SpatialCoordinate(cell)
mu    = Constant(0.1,  "mu")
nu    = Constant(0.01, "nu")
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

mainOp      = galerkinScheme([mainModel==0,DirichletBC(spcU,exact_u,1)])
gradOp      = galerkinOperator(gradModel)
divOp       = galerkinOperator(divModel)
massOp      = galerkinScheme(massModel==0)
preconOp    = galerkinScheme(preconModel==0)

velocity = spcU.interpolate([0,]*spcU.dimRange, name="velocity")
pressure = spcP.interpolate([0], name="pressure")
rhsVelo  = velocity.copy()
rhsPress = pressure.copy()

r      = rhsPress.copy()
d      = rhsPress.copy()
precon = rhsPress.copy()
xi     = rhsVelo.copy()

A = linearOperator(mainOp)
G = linearOperator(gradOp)
D = linearOperator(divOp)
M = linearOperator(massOp)
P = linearOperator(preconOp)

# Note issue with using 'cg' due to
# (a) missing mainOp.setConstraints(velocity)
# (b) no bc for pressure preconder....
solver = {"method":"cg","tolerance":1e-10, "verbose":False}
Ainv   = mainOp.inverseLinearOperator(A,parameters=solver)
Minv   = massOp.inverseLinearOperator(M,parameters=solver)
Pinv   = preconOp.inverseLinearOperator(P,solver)

mainOp(velocity,rhsVelo)
rhsVelo *= -1
G(pressure,xi)
rhsVelo -= xi
mainOp.setConstraints(rhsVelo)
mainOp.setConstraints(velocity)

print("A")
Ainv(rhsVelo, velocity)
D(velocity,rhsPress)
print("B")
Minv(rhsPress, r)

if mainOp.model.nu > 0:
    precon.clear()
    print("C")
    Pinv(rhsPress, precon)
    r *= mainOp.model.mu
    r.axpy(mainOp.model.nu,precon)
d.assign(r)
delta = r.scalarProductDofs(rhsPress)
print("delta:",delta,flush=True)
assert delta >= 0
for m in range(100):
    xi.clear()
    G(d,rhsVelo)
    mainOp.setConstraints(\
        [0,]*grid.dimension, rhsVelo)
    Ainv(rhsVelo, xi)
    D(xi,rhsPress)
    rho = delta /\
       d.scalarProductDofs(rhsPress)
    pressure.axpy(rho,d)
    velocity.axpy(-rho,xi)
    D(velocity, rhsPress)
    Minv(rhsPress,r)
    if mainOp.model.nu > 0:
        precon.clear()
        Pinv(rhsPress,precon)
        r *= mainOp.model.mu
        r.axpy(mainOp.model.nu,precon)
    oldDelta = delta
    delta = r.scalarProductDofs(rhsPress)
    print("delta:",delta,flush=True)
    if delta < 1e-9: break
    gamma = delta/oldDelta
    d *= gamma
    d += r
velocity.plot(colorbar="horizontal")

###############################################################

# Note: the models are virtualize so op.model().nu() will not work on the C++ side
#       so we either need to 'devirtualize' or pass in nu,mu
from dune.generator import algorithm
velocity = spcU.interpolate([0,]*spcU.dimRange, name="velocity")
pressure = spcP.interpolate([0], name="pressure")

# get main forcing and set constraints
mainOp(velocity,rhsVelo)
rhsVelo *= -1
G(pressure,xi)
rhsVelo -= xi
mainOp.setConstraints(rhsVelo)
mainOp.setConstraints(velocity)

uzawa = algorithm.run('uzawa', 'uzawa.hh',
            mainOp.model.nu,mainOp.model.mu,
            mainOp, G, D, Ainv, Minv, Pinv, rhsVelo, xi, rhsPress, r, d,
            precon, velocity, pressure)

fig = pyplot.figure(figsize=(20,10))
velocity.plot(colorbar="horizontal", figure=(fig, 121))
pressure.plot(colorbar="horizontal", figure=(fig, 122))
pyplot.show()
