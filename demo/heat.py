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

ufl2model = duneuflmodel.DuneUFLModel(grid.dimWorld, 1, 'Heat')
u         = ufl2model.trialFunction()
v         = ufl2model.testFunction()
u_n       = ufl2model.coefficient('u_n')

# Crank Nicholson
a = (ufl.inner(u-u_n, v) +
     deltaT * ufl.inner(ufl.grad((u+u_n)/2), ufl.grad(v))) * ufl.dx(0)
ufl2model.generate(a)
heatModel = ufl2model.makeAndImport(grid,name="heat").get()

solution = spc.interpolate(lambda x: [x[0]*(1-x[0])*x[1]*(1-x[1])], name="heat")
grid.writeVTK("heat", pointdata=[solution], number=0)

heatScheme = scheme.create("FemScheme", spc, heatModel, "heat")

u_n = spc.interpolate(solution)
heatModel.setu_n(u_n)

steps = int(1 / deltaT)
for n in range(0,steps):
    # heatModel.setu_n(solution)
    # solution = heatScheme.solve()
    u_n.assign(solution)
    heatScheme.solve( target=solution )
    grid.writeVTK("heat", pointdata=[solution], number=n+1)
