"""Solve the Laplace equation
"""

from mpi4py import MPI

import math
from ufl import *

import dune.ufl
import dune.grid
import dune.fem
import dune.fem.function as gf
import dune.fem.space
import dune.fem.scheme

# dune.fem.create.spaceGenerator.force = True

dune.femmpi.parameter.append("../data/parameter")

grid = dune.grid.create("ALUConform", dune.grid.cartesianDomain([0,0],[1,1],[8,8]), dimgrid=2)
# spc  = dune.fem.create.space("DGONB", grid, dimrange=1, order=2)
spc  = dune.fem.space.create("Lagrange", grid, dimrange=1, order=2)

uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1, field="double")
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )

a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
a = a + 20./(u[0]*u[0]+1.) * v[0] * dx

model = dune.fem.create.ellipticModel(grid, a==0, exact=exact)()

# scheme = dune.fem.create.scheme("DGFemScheme", spc, model,\
scheme = dune.fem.scheme.create("H1", spc, model, "scheme",\
       parameters=\
       {"fem.solver.newton.linabstol": 1e-10,
        "fem.solver.newton.linreduction": 1e-10,
        "fem.solver.newton.verbose": 1,
        "fem.solver.newton.linear.verbose": 1},\
        storage="istl")

exact_gf = grid.function("exact", 5, ufl=exact)
for i in range(2):
    print("solve on level",i)
    uh = scheme.solve()
    def l2error(en,x):
        val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
        return [ val[0]*val[0] ];
    l2error_gf = grid.function( "error", 5, localExpr=l2error )
    error = math.sqrt( grid.l2Norm(l2error_gf) )
    print("size:",grid.size(0),"L2-error:",error)
    grid.writeVTK("laplace", pointdata=[ uh,l2error_gf ])
    grid.hierarchicalGrid.globalRefine(2)
