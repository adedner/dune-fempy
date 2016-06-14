from __future__ import print_function

# from mpi4py import MPI

import dune.fem.grid as grid
import dune.fem.space as space

import math

def testSpace(gridtype, **kwargs):
    grid2d = grid.leafGrid("../data/unitcube-2d.dgf", gridtype, dimgrid=2)
    vtk = grid2d.vtkWriter()
    lagrangespace = space.create("Lagrange", grid2d)

    def expr_global(x):
        return [-(x[1] - 0.5)*math.sin(x[0]*12)]

    gf = grid2d.globalGridFunction("expr_global", expr_global)
    df = lagrangespace.interpolate(gf)

    df.addToVTKWriter(vtk, vtk.PointData)

    df2 = lagrangespace.interpolate([5,3])
    df2.addToVTKWriter(vtk, vtk.CellData)

    df3 = lagrangespace.interpolate(df)
    df3.addToVTKWriter(vtk, vtk.PointData)

    vtk.write("space_demo");

print("ALUGRID")
testSpace("ALUSimplexGrid")
print("YASPGRID A")
testSpace("YaspGrid")
print("YASPGRID B")
testSpace("YaspGrid")
print("END")
