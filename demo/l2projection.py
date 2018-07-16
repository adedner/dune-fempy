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

model = create.model("integrands", grid, a == b)

scheme = create.scheme("galerkin", model, space)
solution, _ = scheme.solve()
grid.writeVTK("l2projection", pointdata=[solution])

scheme.assemble(solution).as_istl.store("l2projection.mm", "matrixmarket")
