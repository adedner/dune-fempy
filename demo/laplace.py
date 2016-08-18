"""Solve the Laplace equation
"""

from mpi4py import MPI

import math
from ufl import *

from dune.models.elliptic import importModel
import dune.ufl
import dune.fem

dune.femmpi.parameter.append("../data/parameter")

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

f = cos(2*math.pi*x[0])*cos(2*math.pi*x[1])

a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
b = f * v[0] * dx
model = importModel(grid, a == b).get()

scheme = dune.fem.create.scheme("FemScheme", spc, model, "scheme",\
   {"fem.solver.newton.linabstol": 1e-9,
    "fem.solver.newton.linreduction": 1e-9,
    "fem.solver.newton.verbose": 1,
    "fem.solver.newton.linear.verbose": 1})
grid.writeVTK("laplace", pointdata=[scheme.solve()])
# print(str(dune.femmpi.parameter))
print(model.dimRange())
