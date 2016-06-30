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

a = (inner(grad(u), grad(v)) + inner(u,v)) * dx(0)
b = sin(2*math.pi*x[0])*sin(2*math.pi*x[1]) * v[0] * dx(0)

model = dune.models.elliptic.compileUFL(a == b, dimRange=1)

scheme = dune.fem.create.scheme("FemScheme", spc, dune.models.elliptic.importModel("laplacemodel", grid, model).get(), "scheme")
grid.writeVTK("laplace", pointdata=[scheme.solve()])
