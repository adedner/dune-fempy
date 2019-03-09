from __future__ import print_function

import math

import ufl
import dune.ufl
import dune.fem
import dune.fem.function as gf

import dune.create as create

dune.fem.parameter.append("../data/parameter")

def getUflSpace(grid,spc):
    # test different versions of constructing a uflSpace
    uflSpace = dune.ufl.Space(grid.dimGrid, 1, field="double")
    x = ufl.SpatialCoordinate(uflSpace.cell())
    uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1, field="double")
    x = ufl.SpatialCoordinate(uflSpace.cell())
    uflSpace = dune.ufl.Space(grid, 1, field="double")
    x = ufl.SpatialCoordinate(uflSpace.cell())
    uflSpace = dune.ufl.Space(spc)
    x = ufl.SpatialCoordinate(uflSpace.cell())
    return uflSpace, x

def compute():
    grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8]), dimgrid=2)
    spc  = dune.create.space("Lagrange", grid, dimRange=1, order=2, storage="istl")

    uflSpace, x = getUflSpace(grid,spc)

    exact = ufl.as_vector( [ufl.cos(2.*ufl.pi*x[0])*ufl.cos(2.*ufl.pi*x[1])] )
    exact_gf = create.function("ufl", grid, "exact", 5, exact)

    uh = spc.interpolate( exact_gf, name="solution")

    # Approach 1
    def l2error(en,x):
        val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
        return [ val[0]*val[0] ];
    l2error_gf = create.function("local", grid, "error", 5, l2error )
    error = math.sqrt( l2error_gf.integrate()[0] )
    print("Approach 1:",error)

    # Approach 2
    uh_coeff = ufl.Coefficient(uflSpace)
    l2error_gf = create.function("ufl", grid, "error", 5,
            ufl.as_vector([(exact[0]-uh_coeff[0])**2]), coefficients={uh_coeff:uh} )
    error = math.sqrt( l2error_gf.integrate()[0] )
    print("Approach 2:",error)

    # Approach 3
    l2error_gf = create.function("ufl", grid, "error", 5, (exact-uh)**2 )
    error = math.sqrt( l2error_gf.integrate()[0] )
    print("Approach 4:",error)

compute()
