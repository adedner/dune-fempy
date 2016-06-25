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

deltaT = 0.01

grid = fem.leafGrid(dgf, "ALUSimplexGrid", dimgrid=2)
spc = space.create( "Lagrange", grid, dimrange=1, polorder=2)

# why dimWorld?
model = duneuflmodel.DuneUFLModel(grid.dimWorld, 1, 'Heat')
u = model.trialFunction()
v = model.testFunction()
u_n = model.coefficient('u_n')

a = (ufl.inner(u, v) - deltaT * ufl.inner(ufl.grad(u), ufl.grad(v))) * ufl.dx(0)
b = ufl.inner(u_n, v)
model.generate(a,b)
heatModel = model.makeAndImport(grid,name="heat").get()

solution = spc.interpolate(lambda x: [x[0]*(1-x[0])*x[1]*(1-x[1])], name="solution")
grid.writeVTK("heat", pointdata=[solution], number=0)

heatScheme = scheme.create("FemScheme", spc, heatModel, "laplace")

steps = int(1 / deltaT)
for n in range(0,steps):
    heatScheme.setu_n(solution)
    solution = heatScheme.solve()
    grid.writeVTK("heat", pointdata=[solution], number=n+1)
