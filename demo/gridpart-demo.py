from __future__ import print_function
import math
from mpi4py import MPI

import dune.fem as fem
import dune.fem.gridpart as gridpart

def testGridPart(gridtype):
    grid2d = fem.leafGrid("../data/unitcube-2d.dgf", gridtype, dimgrid=2)

    def expr_global(x):
        return [x[0]*(x[0]+1),(x[0]+1.)*x[1]] # ,math.sin(x[0]*x[1]*2*math.pi)] # problem in vtk with dimensionworld increase in geogp

    gf = grid2d.globalGridFunction("expr_global", expr_global)
    df = grid2d.interpolate(gf, space="Lagrange", name="test")

    geogp = gridpart.create("Geometry", df )
    vtk = geogp.vtkWriter()
    gfnew = geogp.globalGridFunction("global", expr_global)

    gfnew.addToVTKWriter(vtk, vtk.PointData)
    vtk.write("gridpart_demo");

print("YASPGRID B")
testGridPart("YaspGrid")
print("END")
