use_codegen = True

import time
from dune.grid import structuredGrid
from dune.fem import parameter
import dune.create as create
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle, exp,\
                dx, grad, inner, as_vector, replace, sqrt
from dune.ufl import NamedConstant

parameter.append({"fem.verboserank": -1})

grid = structuredGrid([0, 0], [1, 1], [40, 40])

order     = 4
dimR      = 3
quadOrder = 2*order+3
if use_codegen:
    space = create.space('lagrange', grid, dimrange=dimR, order=order,
                         interiorQuadratureOrders=[quadOrder],\
                         skeletonQuadratureOrders=[quadOrder] )
else:
    space = create.space('lagrange', grid, dimrange=dimR, order=order)

x = SpatialCoordinate(space)

initial = 1/2*(x[0]**2 + x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1

u_h = space.interpolate(as_vector([initial,]*dimR), name='u_h')
u_h_n = u_h.copy(name="previous")

u = TrialFunction(space)
v = TestFunction(space)
dt = NamedConstant(triangle, "dt")    # time step
t  = NamedConstant(triangle, "t")     # current time

abs_du = sqrt(inner(grad(u), grad(u)))
K = 2/(1 + sqrt(1 + 4*abs_du))
a = (inner((u - u_h_n)/dt, v) + inner(K*grad(u), grad(v)))*dx
exact = as_vector( [exp(-2*t)*(initial - 1) + 1,]*dimR )
b = replace(a, {u: exact})

solverParam = {"fem.solver.newton.verbose": 0,
               "fem.solver.newton.linear.verbose": 0}
scheme = create.scheme("galerkin", a==b, space, solver='cg',
                       parameters = solverParam)
scheme.setQuadratureOrders(quadOrder,quadOrder)

scheme.model.dt = 0.05

grid.writeVTK('initial', pointdata={'initial': initial})

start = time.time()
t = 0
while t < 1.0:
    scheme.model.t = t
    u_h_n.assign(u_h)
    # scheme.solve(target=u_h)
    scheme.assemble(u_h_n)
    t += scheme.model.dt
print("time loop:",time.time()-start)

grid.writeVTK('forchheimer', pointdata=[u_h])
