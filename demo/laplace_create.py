"""Solve the Laplace equation
"""
from __future__ import print_function

import math
from ufl import *

import dune.ufl
import dune.fem
import dune.fem.function as gf

import dune.create as create

dune.fem.parameter.append("../data/parameter")

def compute():
    uflSpace = dune.ufl.Space((2,2), 1, field="double")
    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    x = SpatialCoordinate(uflSpace.cell())

    exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )

    a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
    a = a + 20./(u[0]*u[0]+1.) * v[0] * dx

    eqn = a==0
    domain = dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8])
    parameters = {"fem.solver.newton.linabstol": 1e-10,
                "fem.solver.newton.linreduction": 1e-10,
                "fem.solver.newton.verbose": 0,
                "fem.solver.gmres.restart": 50,
                "fem.solver.newton.linear.verbose": 0}
    # using both grid and view (AdaptivLeaf<ALU>)
    scheme = create.scheme("h1", storage="istl",
               model="elliptic", equation=eqn, exact=exact, dirichlet={ 1:exact },
               space="Lagrange",order=2,
               grid="ALUCube", constructor=domain,dimgrid=2,
               view="adaptive",
               parameters=parameters)
    # using only grid (LeafGrid<ALU>)
    scheme = create.scheme("h1", storage="istl",
               model="elliptic", equation=eqn, exact=exact, dirichlet={ 1:exact },
               space="Lagrange",order=2,
               grid="ALUCube", constructor=domain,dimgrid=2,
               parameters=parameters)
    # using only view (LeafGrid<ALU>)
    scheme = create.scheme("h1", storage="istl",
               model="elliptic", equation=eqn, exact=exact, dirichlet={ 1:exact },
               space="Lagrange",order=2,
               view="ALUCube", constructor=domain,dimgrid=2,
               parameters=parameters)

    grid = scheme.space.grid
    # uh = create.discretefunction("istl",scheme.space,"solution")

    exact_gf = create.function("ufl", grid, "exact", 5, exact)
    for i in range(2):
        print("solve on level",i, "number of dofs=",grid.size(2))
        uh = scheme.solve()
        def l2error(en,x):
            val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
            return [ val[0]*val[0] ];
        l2error_gf = create.function("local", grid, "error", 5, l2error )
        error = math.sqrt( l2error_gf.integrate()[0] )
        print("size:",grid.size(0),"L2-error:",error)
        grid.writeVTK("laplace", pointdata=[ uh,l2error_gf ])
        grid.hierarchicalGrid.globalRefine(1)
        print("  ***** REFINE GRID ***** ")

compute()
