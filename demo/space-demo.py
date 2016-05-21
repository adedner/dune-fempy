from __future__ import print_function

# from mpi4py import MPI

import dune.fem.grid as grid
import dune.fem.space as space

import math

yaspgrid = grid.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)
lagrangespace = space.create("Lagrange", yaspgrid)

def expr_global(x):
    return [-(x[1] - 0.5)*math.sin(x[0]*12)]

gf = yaspgrid.globalGridFunction("expr_global", expr_global)
df = lagrangespace.interpolate(gf)

vtk_yaspgrid = yaspgrid.vtkWriter()
df.addToVTKWriter(vtk_yaspgrid, vtk_yaspgrid.DataType.PointData)
df2 = lagrangespace.interpolate([5])
df2.addToVTKWriter(vtk_yaspgrid, vtk_yaspgrid.DataType.CellData)
vtk_yaspgrid.write("space_demo");
