from __future__ import print_function

import math

import dune.common as common
from dune.fem.gridpart import create as gridPart
#from dune.fem.gridpart.geometry import create as geometryGridPart
#from dune.fem.gridpart.filtered import create as filteredGridPart

import dune.create as create

def testGeometryGridPart(grid, prefix):
    t = 0
    def expr_global(x):
        return [x[0]*(x[0]+1),(x[0]+1.)*x[1]*math.sin(0.1+2.*math.pi*t)] # ,math.sin(x[0]*x[1]*2*math.pi)] # problem in vtk with dimensionworld increase in geogp

    gf = grid.function("expr_global", order=1, globalExpr=expr_global)
    df = grid.interpolate(gf, space="Lagrange", name="test")

    geogp = gridPart("geometry", df)
    #geogp = geometryGridPart(df)
    vtk = geogp.vtkWriter()
    gfnew = geogp.function("global", order=1, globalExpr=expr_global)
    gfnew.addToVTKWriter(vtk, common.DataType.PointData)

    dt = 0.01
    count = 0
    while t<1:
        t += dt
        count += 1
        df.interpolate(gf)
        vtk.write(prefix + str(count));

def testGridPart(gridtype):
    grid = create.grid(gridtype, "../data/unitcube-2d.dgf", dimgrid=2)
    testGeometryGridPart(grid, "gridpart_demo")

    #subGrid = filteredGridPart(grid, lambda e: (e.geometry.center - [0.5, 0.5]).two_norm < 0.25)
    # subGrid = gridPart("filtered", grid, lambda e: (e.geometry.center - [0.5, 0.5]).two_norm < 0.25)
    # testGeometryGridPart(subGrid, "gridpart_demo_sub")

print("YASPGRID B")
testGridPart("Yasp")
print("END")
