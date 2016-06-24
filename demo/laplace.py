from __future__ import print_function
import math
from mpi4py import MPI

import ufl
import dune.models.femufl as duneuflmodel
import dune.fem as fem
import dune.fem.space as space
import dune.fem.scheme as scheme

dgf = """DGF

INTERVAL
0  0
1  1
16 16
#
"""

grid = fem.leafGrid(dgf, "ALUSimplexGrid", dimgrid=2)
spc = space.create( "Lagrange", grid, dimrange=1, polorder=2)

# why dimWorld?
model = duneuflmodel.DuneUFLModel(grid.dimWorld, 1)
u = model.trialFunction()
v = model.testFunction()
x = model.spatialCoordinate()

a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx(0)
b = ufl.sin(x[0])*ufl.sin(x[1]) * v[0] * ufl.dx(0)
model.generate(a,b)
laplaceModel = model.makeAndImport(grid,name="laplace").get()

laplace = scheme.create("FemScheme", spc, laplaceModel, "laplace")

grid.writeVTK("laplace", pointdata=[laplace.solve()])
