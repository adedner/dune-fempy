from __future__ import print_function
import math,sympy
from mpi4py import MPI

from ufl import *
import dune.models.femufl as duneuflmodel
import dune.fem as fem

grid = fem.leafGrid(fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc = fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

# why dimWorld?
ufl2model = duneuflmodel.DuneUFLModel(grid.dimWorld, 1)
u = ufl2model.trialFunction()
v = ufl2model.testFunction()
x = ufl2model.spatialCoordinate()

#f = sin(2*math.pi*x[0])*sin(2*math.pi*x[1])
f = (x[0] + x[1])

a = (inner(grad(u), grad(v)) + inner(u,v)) * dx(0)
b = f * v[0] * dx(0)

ufl2model.generate(a-b)
model = ufl2model.makeAndImport(grid,name="laplace").get()

laplace = fem.create.scheme("FemScheme", spc, model, "laplace")
grid.writeVTK("laplace_1", pointdata=[laplace.solve()])

#ufl2model.clear()
#
#bnd   = 100.*ufl2model.x0*ufl2model.x1*(1-ufl2model.x0)*(1-ufl2model.x1)
#exact = [sympy.atan(bnd*bnd)]
#ufl2model.generate(a,exact=exact)
#model = ufl2model.makeAndImport(grid,name="laplace",exact=exact).get()
#
#laplace = scheme.create("FemScheme", spc, model, "laplace")
#grid.writeVTK("laplace_2", pointdata=[laplace.solve()])
