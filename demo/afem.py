from __future__ import print_function
import math,sympy
from mpi4py import MPI

import ufl
import dune.models.femufl as duneuflmodel
import dune.fem as fem
import dune.fem.space as space
import dune.fem.scheme as scheme
import dune.fem.function as function

dgf = """DGF

INTERVAL
-1 -1
1  1
16 16
#
"""

grid = fem.leafGrid("../data/unitcube-2d.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc = space.create( "Lagrange", grid, dimrange=1, polorder=2)

# why dimWorld?
model     = duneuflmodel.DuneUFLModel(2,1) # grid.dimWorld, 1)
u         = model.trialFunction()
v         = model.testFunction()
x         = model.spatialCoordinate()

a = ( ufl.inner(ufl.grad(u), ufl.grad(v)) + ufl.inner(u,v) ) * ufl.dx(0)
# L = (16.-10.*ufl.atan( 1000.*(x[0]*x[0]+x[1]*x[1])))*v[0] * ufl.dx(0)
# func  = 16.-10.*sympy.atan(1000.*(model.x0*model.x1+model.x1*model.x1))
bnd   = 100.*model.x0*model.x1*(1-model.x0)*(1-model.x1)
exact = [sympy.atan(bnd*bnd)]
model.generate(a,exact=exact)
laplaceModel = model.makeAndImport(grid,name="laplace",exact=exact).get()

laplace = scheme.create("FemScheme", spc, laplaceModel, "laplace")
uh = laplace.solve()
grid.writeVTK("laplace", pointdata=[uh])

count = 1
tol = 0.3
while count<10:
    [estimate,marked] = laplace.mark(uh,tol)
    if marked == False:
        break
    grid.hierarchicalGrid.adapt([uh])
    grid.hierarchicalGrid.loadBalance([uh])
    # grid.hierarchicalGrid.globalRefine(2)
    laplace.solve( target=uh )
    grid.writeVTK("laplace",
            pointdata=[uh],celldata=[grid.localGridFunction("level", function.Levels())], number=str(count) )
    count += 1
    print(count, ": size=",grid.size(0), "estimate=",estimate,"error=",laplace.error(uh))
