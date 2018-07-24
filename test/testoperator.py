import time
from dune.grid import structuredGrid
from dune.fem import parameter
from dune.fem.function import integrate
import dune.create as create
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle, exp,\
                dx, grad, inner, as_vector, replace, sqrt
from dune.ufl import NamedConstant

parameter.append({"fem.verboserank": -1})

grid = structuredGrid([0, 0], [1, 1], [40, 40])

order     = 1
dimR      = 1
quadOrder = 2*order+3
spaceName = "lagrange"
space = create.space(spaceName, grid, dimrange=dimR, order=order)

arg   = space.interpolate(as_vector([1,]*dimR), name='arg')
destA = space.interpolate(as_vector([0,]*dimR), name='destA')
destB = space.interpolate(as_vector([0,]*dimR), name='destB')

u  = TrialFunction(space)
v  = TestFunction(space)
a  = ( inner(u, v) + inner(grad(u), grad(v)) ) * dx
a  = inner(u, v) * dx
op = create.operator("galerkin", a, space)

A = op.assemble(arg)
op(arg,destA)
A(arg,destB)
# grid.writeVTK('optest', pointdata=[destA,destB])
err = integrate(grid, (destA-destB)**2, 5)
print("error=",err)
