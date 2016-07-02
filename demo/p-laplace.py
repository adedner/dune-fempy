from __future__ import print_function

from mpi4py import MPI

import math
from ufl import *

import dune.models.elliptic
import dune.ufl
import dune.fem

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

uflSpace = dune.ufl.Space(grid.dimGrid, 1, grid.dimWorld)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

d = 0.001
p = 1.7

rhs = (x[0] + x[1]) * v[0]

a = (pow(d + inner(grad(u), grad(u)), (p-2)/2)*inner(grad(u), grad(v)) + inner(u, v)) * dx(0) + 10*inner(u, v) * ds(0)
#b = sin(2*math.pi*x[0])*sin(2*math.pi*x[1]) * v[0] * dx(0)
b = rhs * dx(0) + 10*rhs * ds(0)

model = dune.models.elliptic.compileUFL(a == b)

scheme = dune.fem.create.scheme("FemScheme", spc, dune.models.elliptic.importModel(grid, model).get(), "scheme")
grid.writeVTK("p-laplace", pointdata=[scheme.solve()])
