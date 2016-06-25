from __future__ import print_function
import math,sympy
from mpi4py import MPI

import ufl
import dune.models.femufl as duneuflmodel
import dune.fem as fem

dgf = """DGF

INTERVAL
 0  0
 1  1
 5  5
#
"""

grid = fem.leafGrid(dgf, "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc  = fem.create.space( "Lagrange", grid, dimrange=1, polorder=2)

ufl2model = duneuflmodel.DuneUFLModel(grid.dimWorld, 1)
u         = ufl2model.trialFunction()
v         = ufl2model.testFunction()
x         = ufl2model.spatialCoordinate()
x0        = ufl2model.x0
x1        = ufl2model.x1

a = ( ufl.inner(ufl.grad(u), ufl.grad(v)) + ufl.inner(u,v) ) * ufl.dx(0)
bnd   = 100.*x0*x1*(1-x0)*(1-x1)
exact = [sympy.atan(bnd*bnd)]
ufl2model.generate(a,exact=exact)
laplaceModel = ufl2model.makeAndImport(grid,name="laplace",exact=exact).get()

laplace = fem.create.scheme("FemScheme", spc, laplaceModel, "afem")
uh = laplace.solve()

count = 0
tol = 0.3
while count < 20:
    [estimate, marked] = laplace.mark(uh, tol)
    grid.writeVTK("afem", pointdata=[uh], celldata=[grid.levelFunction()], number=str(count) )
    print(count, ": size=",grid.size(0), "estimate=",estimate,"error=",laplace.error(uh))
    if marked == False or estimate < tol:
        break
    grid.hierarchicalGrid.adapt([uh])
    grid.hierarchicalGrid.loadBalance([uh])
    laplace.solve( target=uh )
    count += 1
