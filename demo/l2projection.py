from __future__ import print_function

import math
from ufl import *

from dune.grid import cartesianDomain
from dune.alugrid import aluConformGrid

from dune.ufl import Space as UFLSpace

import dune.create as create

domain = cartesianDomain([0, 0], [1, 1], [8, 8])
grid = aluConformGrid(domain, dimgrid=2)

lagrangeSpace = create.space("Lagrange", grid, dimrange=1, order=2, storage="istl")
dgSpace = create.space("DGONB", grid, dimrange=1, order=2, storage="istl")
space = create.space("combined", lagrangeSpace, dgSpace)

uflSpace = UFLSpace(2, 2)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

a = inner(u, v) * dx
f = sin(math.pi*x[0]) * sin(math.pi * x[1])
b = inner(as_vector([f, f]), v) * dx

scheme = create.scheme("h1", a==b, space)
solution = space.interpolate([0],name="solution")
scheme.solve(target=solution)
grid.writeVTK("l2projection", pointdata=[solution])

from dune.fem.operator import linear
op = linear(space)
scheme.jacobian(solution,op)
op.as_istl.store("l2projection.mm", "matrixmarket")
