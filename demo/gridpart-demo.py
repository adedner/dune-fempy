from __future__ import print_function
import math
from mpi4py import MPI

import dune.fem as fem
import dune.fem.gridpart as gridpart

def testGridPart(gridtype):
    grid2d = fem.leafGrid("../data/unitcube-2d.dgf", gridtype, dimgrid=2)

    t = 0
    def expr_global(x):
        return [x[0]*(x[0]+1),(x[0]+1.)*x[1]*math.sin(0.1+2.*math.pi*t)] # ,math.sin(x[0]*x[1]*2*math.pi)] # problem in vtk with dimensionworld increase in geogp

    gf = grid2d.globalGridFunction("expr_global", expr_global)
    df = grid2d.interpolate(gf, space="Lagrange", name="test")

    geogp = gridpart.create("Geometry", df )
    vtk = geogp.vtkWriter()
    gfnew = geogp.globalGridFunction("global", expr_global)
    gfnew.addToVTKWriter(vtk, vtk.PointData)

    dt = 0.01
    count = 0
    while t<1:
        t += dt
        count += 1
        df.interpolate(gf)
        vtk.write("gridpart_demo"+str(count));

print("YASPGRID B")
testGridPart("YaspGrid")
print("END")
