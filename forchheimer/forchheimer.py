import time
import numpy
from dune.grid import structuredGrid

from dune.fem import parameter
parameter.append({"fem.verboserank": 0})

grid = structuredGrid([0, 0], [1, 1], [4, 4])

import dune.create as create
space = create.space('lagrange', grid, dimrange=1, order=2)

from ufl import SpatialCoordinate
x = SpatialCoordinate(space)

initial = 1/2*(x[0]**2 + x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1

u_h = space.interpolate(initial, name='u_h')
u_h_n = u_h.copy(name="previous")

from ufl import TestFunction, TrialFunction, triangle, exp,\
                dx, grad, inner, as_vector, replace, sqrt
from dune.ufl import NamedConstant
u = TrialFunction(space)
v = TestFunction(space)
dt = NamedConstant(triangle, "dt")    # time step
t  = NamedConstant(triangle, "t")     # current time

abs_du = sqrt(inner(grad(u), grad(u)))
K = 2/(1 + sqrt(1 + 4*abs_du))
a = (inner((u - u_h_n)/dt, v) + inner(K*grad(u), grad(v)))*dx
exact = as_vector( [exp(-2*t)*(initial - 1) + 1] )
b = replace(a, {u: exact})

solverParam = {"fem.solver.newton.verbose": 1,
               "fem.solver.newton.linear.verbose": 1}
scheme = create.scheme("h1", space, a == b, solver='cg', parameters = solverParam)

scheme.model.dt = 0.05

grid.writeVTK('initial', pointdata={'initial': initial})

time = 0
while time < 1.0:
    scheme.model.t = time
    u_h_n.assign(u_h)
    scheme.solve(target=u_h)
    time += scheme.model.dt

grid.writeVTK('forchheimer', pointdata=[u_h])
