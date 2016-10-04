"""Solve the Laplace equation
"""
from __future__ import print_function

import math
from ufl import *

import dune.ufl
import dune.fem
import dune.fem.function as gf

from dune.fem.view import adaptiveLeafGridView

import dune.create as create

dune.fem.parameter.append("../data/parameter")

def compute():
    # grid = create.grid("SPIsotropic", dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8]), dimgrid=2)
    grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8]), dimgrid=2)
    grid = adaptiveLeafGridView(grid)

    # spc  = dune.fem.create.space("DGONB", grid, dimrange=1, order=2)
    spc  = dune.create.space("Lagrange", grid, dimrange=1, order=2)

    uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1, field="double")
    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    x = SpatialCoordinate(uflSpace.cell())

    exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )

    a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
    a = a + 20./(u[0]*u[0]+1.) * v[0] * dx

    model = create.model("elliptic", grid, a==0, exact=exact, dirichlet={ 1:exact } )

    # scheme = create.scheme("DGFemScheme", spc, model,\
    scheme = create.scheme("h1", spc, model, "scheme",\
           parameters=\
           {"fem.solver.newton.linabstol": 1e-10,
            "fem.solver.newton.linreduction": 1e-10,
            "fem.solver.newton.verbose": 0,
            "fem.solver.newton.linear.verbose": 0},\
            storage="istl")

    exact_gf = create.function("ufl", grid, "exact", 5, exact)
    for i in range(2):
        print("solve on level",i)
        uh = scheme.solve()
        def l2error(en,x):
            val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
            return [ val[0]*val[0] ];
        l2error_gf = create.function("local", grid, "error", 5, l2error )
        error = math.sqrt( l2error_gf.integrate()[0] )
        print("size:",grid.size(0),"L2-error:",error)
        grid.writeVTK("laplace", pointdata=[ uh,l2error_gf ])
        grid.hierarchicalGrid.globalRefine(2)

compute()
