from __future__ import print_function

import math

import dune.fem.space as space
# import dune.fem.discretefunction as discfunc

import dune.create as create

def testSpace(grid2d, spacetype, dimRange, order):
    lagrangespace = space.create(spacetype, grid2d, dimrange=dimRange, order=order)

    # would work but perhaps not desired
    # df0 = discfunc.create("Adaptive",lagrangespace,name="test")
    # df0.clear();
    # df00 = discfunc.create("Adaptive",lagrangespace,name="testA")
    # df00.clear();
    # df0.addToVTKWriter(vtk,vtk.PointData)
    # df00.addToVTKWriter(vtk,vtk.PointData)

    expr_global = lambda x: [-(x[1] - 0.5)*math.sin(x[0]*12)]*dimRange

    gf = grid2d.function("expr_global", order=1, globalExpr=expr_global)
    df1 = lagrangespace.interpolate( gf, name="test" )
    df2 = lagrangespace.interpolate( [5]*dimRange )
    df3 = lagrangespace.interpolate( df1, name="copy", storage="istl" )
    df4 = lagrangespace.interpolate( lambda x:
            [(x-[0.5,0.5]).infinity_norm,]*dimRange, name="radius") # , storage="Numpy"  )
    lagrangespace=0
    grid2d.writeVTK("space_demo", pointdata=[gf,df1,df2,df3,df4])

def test(gridtype):
    grid2d = create.grid(gridtype, "../data/unitcube-2d.dgf", dimgrid=2)
    print("Lagrange(1,1)")
    testSpace(grid2d, "Lagrange",1,1)
    print("Lagrange(1,2)")
    testSpace(grid2d, "Lagrange",1,2)
    print("DGONB(2,0)")
    testSpace(grid2d, "DGONB",2,0)
    print("DGONB(1,2)")
    testSpace(grid2d, "DGONB",1,2)
    print("P1Bubble(2)")
    testSpace(grid2d, "P1Bubble",2,1)


print("ALUSimplex")
test("ALUSimplex")
print("YASPGRID A")
test("Yasp")
print("YASPGRID B")
test("Yasp")
print("END")
