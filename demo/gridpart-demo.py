from __future__ import print_function

import math

import dune.common as common
from dune.fem.gridpart import create as gridPart
#from dune.fem.gridpart.geometry import create as geometryGridPart
#from dune.fem.gridpart.filtered import create as filteredGridPart
from dune.fem.view import geometryGridView, filteredGridView

import dune.create as create

def testGeometryGridView(grid, prefix):
    t = 0
    def expr_global(x):
        return [x[0]*(x[0]+1),(x[0]+1.)*x[1]*math.sin(0.1+2.*math.pi*t)] # ,math.sin(x[0]*x[1]*2*math.pi)] # problem in vtk with dimensionworld increase in geogp

    gf = create.function("global", grid, "coordinates", 1, expr_global)
    spc = create.space("Lagrange", grid, dimrange=2, order=1)
    df = spc.interpolate(gf, name="test")

    geogrid = geometryGridView(df)
    gfnew = create.function("global", geogrid, "expression", 1, expr_global)

    dt = 0.01
    count = 0
    while t < 1:
        t += dt
        count += 1
        df.interpolate(gf)
        geogrid.writeVTK(prefix, pointdata=[gfnew], number=count)

def testGridView(gridtype):
    grid = create.grid(gridtype, "../data/unitcube-2d.dgf", dimgrid=2)
    testGeometryGridView(grid, "gridpart-demo")

    subGrid = filteredGridView(grid, lambda e: (e.geometry.center - [0.5, 0.5]).two_norm < 0.25)
    testGeometryGridView(subGrid, "gridpart-demo-sub")

testGridView("Yasp")
