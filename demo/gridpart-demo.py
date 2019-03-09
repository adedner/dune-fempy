from __future__ import print_function

import math

from dune.fem.view import geometryGridView, filteredGridView

import dune.create as create

from dune.plotting import block, disable

def plot(grid):
    if disable: return
    try:
        from matplotlib import pyplot
        from numpy import amin, amax, linspace

        triangulation = grid.triangulation(4)

        pyplot.gca().set_aspect('equal')
        pyplot.triplot(grid.triangulation(), antialiased=True, linewidth=0.2, color='black')
        pyplot.show(block=block)
    except ImportError:
        pass

def testGeometryGridView(grid, prefix):
    t = 0
    def expr_global(x):
        return [x[0]*(x[0]+1),(x[0]+1.)*x[1]*math.sin(0.1+2.*math.pi*t)] # ,math.sin(x[0]*x[1]*2*math.pi)] # problem in vtk with dimensionworld increase in geogp

    gf = create.function("global", grid, "coordinates", 1, expr_global)
    spc = create.space("lagrange", grid, dimRange=2, order=1)
    df = spc.interpolate(gf, name="test")

    geogrid = geometryGridView(df)
    gfnew = create.function("global", geogrid, "expression", 1, expr_global)

    vtk = geogrid.sequencedVTK(prefix, pointdata=[gfnew])
    vtk()

    dt = 0.01
    while t < 1:
        t += dt
        df.interpolate(gf)
        vtk()
    plot(geogrid)

def testGridView(gridtype):
    grid = create.grid(gridtype, "../data/unitcube-2d.dgf", dimgrid=2)
    testGeometryGridView(grid, "gridpart-demo")

    subGrid = filteredGridView(grid, lambda e: (e.geometry.center - [0.5, 0.5]).two_norm < 0.25, 1)
    testGeometryGridView(subGrid, "gridpart-demo-sub")
    plot(subGrid)

testGridView("Yasp")
