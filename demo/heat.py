from __future__ import print_function
import math
from mpi4py import MPI
import dune.fem as fem

import ufl
import dune.models.femufl as duneuflmodel

dgf = """DGF

INTERVAL
0  0
1  1
16 16
#
"""

deltaT = 0.01

grid = fem.leafGrid(dgf, "ALUSimplexGrid", dimgrid=2)
spc = fem.create.space( "Lagrange", grid, dimrange=1, polorder=2)

ufl2model = duneuflmodel.DuneUFLModel(grid.dimWorld, 1, 'Heat')
u         = ufl2model.trialFunction()
v         = ufl2model.testFunction()
u_n       = ufl2model.coefficient('u_n')

# Crank Nicholson
theta = 0.5
a = (ufl.inner(u-u_n, v) +
     deltaT * ufl.inner(ufl.grad(theta*u+(1-theta)*u_n), ufl.grad(v))) * ufl.dx(0)
ufl2model.generate(a)
heatModel = ufl2model.makeAndImport(grid,name="heat").get()

solution = spc.interpolate(lambda x: [ math.atan( (10.*x[0]*(1-x[0])*x[1]*(1-x[1]))**2 ) ], name="heat")
heatScheme = fem.create.scheme("FemScheme", spc, heatModel, "heat")

grid.writeVTK("heat", pointdata=[solution], number=0)
# vtk = grid.writeVTK("heat", pointdata=[solution], number=0)
# u_n = spc.interpolate(solution)
# heatModel.setu_n(u_n)

steps = int(1 / deltaT)
for n in range(1,steps+1):
    heatModel.setu_n(solution)
    solution = heatScheme.solve()
    grid.writeVTK("heat", pointdata=[solution], number=n)
    # u_n.assign(solution)
    # heatScheme.solve( target=solution )
    # vtk.write("heat", n)
