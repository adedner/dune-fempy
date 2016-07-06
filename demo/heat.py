from __future__ import print_function

from mpi4py import MPI

import math
from ufl import *

import dune.models.elliptic
import dune.ufl
import dune.fem

# Crank Nicholson
theta = 0.5
deltaT = 0.01

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
u_n = Coefficient(uflSpace)

a = (inner(u - u_n, v) + deltaT * inner(grad(theta*u + (1-theta)*u_n), grad(v))) * dx(0)

model = dune.models.elliptic.importModel(grid, dune.models.elliptic.compileUFL(a == 0)).get()
scheme = dune.fem.create.scheme("FemScheme", spc, model, "scheme")

solution = spc.interpolate(lambda x: [ math.atan( (10.*x[0]*(1-x[0])*x[1]*(1-x[1]))**2 ) ], name="u")
grid.writeVTK("heat", pointdata=[solution], number=0)

old_solution = spc.interpolate(lambda x: [ math.atan( (10.*x[0]*(1-x[0])*x[1]*(1-x[1]))**2 ) ], name="u_n")
model.setCoefficient(u_n.count(), old_solution)

steps = int(1 / deltaT)
for n in range(1,steps+1):
    old_solution.assign(solution)
    scheme.solve(target=solution)
    grid.writeVTK("heat", pointdata=[solution], number=n)
