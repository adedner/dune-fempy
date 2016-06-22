from __future__ import print_function

from mpi4py import MPI

import dune.fem.grid as grid
import dune.fem.space as space
# import dune.fem.discretefunction as discfunc

import math

def testSpace(gridtype, dimRange):
    grid2d = grid.leafGrid("../data/unitcube-2d.dgf", gridtype, dimgrid=2)
    vtk = grid2d.vtkWriter()
    lagrangespace = space.create("Lagrange", grid2d, dimrange=dimRange)

    # would work but perhaps not desired
    # df0 = discfunc.create("Adaptive",lagrangespace,name="test")
    # df0.clear();
    # df00 = discfunc.create("Adaptive",lagrangespace,name="testA")
    # df00.clear();
    # df0.addToVTKWriter(vtk,vtk.PointData)
    # df00.addToVTKWriter(vtk,vtk.PointData)

    def expr_global1(x):
        return [-(x[1] - 0.5)*math.sin(x[0]*12)]
    def expr_global2(x):
        return [-(x[1] - 0.5)*math.sin(x[0]*12),x[0]*x[1]]

    if dimRange==1:
        gf = grid2d.globalGridFunction("expr_global", expr_global1)
    else:
        gf = grid2d.globalGridFunction("expr_global", expr_global2)
    df1 = lagrangespace.interpolate( gf, name="test" )
    df2 = lagrangespace.interpolate( [5,3] ) # , storage="Numpy" ) # , name="53" )
    df3 = lagrangespace.interpolate( df1, name="copy", storage="Istl" )
    df4 = lagrangespace.interpolate( lambda x: [(x-[0.5,0.5]).infinity_norm,]*dimRange, name="radius"  )
    lagrangespace=0

    gf.addToVTKWriter(vtk, vtk.PointData)
    df1.addToVTKWriter(vtk, vtk.PointData)
    df2.addToVTKWriter(vtk, vtk.CellData)
    df3.addToVTKWriter(vtk, vtk.PointData)
    df4.addToVTKWriter(vtk, vtk.PointData)

    gf=0
    df1=0
    df2=0
    df3=0
    df4=0

    vtk.write("space_demo");

print("ALUGRID")
testSpace("ALUSimplexGrid",1)
print("YASPGRID A")
testSpace("YaspGrid",1)
print("YASPGRID B")
testSpace("YaspGrid",2)
print("END")
