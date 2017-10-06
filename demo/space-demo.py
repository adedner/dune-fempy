from __future__ import print_function

import math

print("import dune.fem.space")
import dune.fem.space as space
# import dune.fem.discretefunction as discfunc

print("import dune.create")
import dune.create as create

def testSpace(grid2d, spacetype, dimRange, order):
    print("HALLO A")
    lagrangespace = create.space(spacetype, grid2d, dimrange=dimRange, order=order, storage="istl")
    print("HALLO B")

    # would work but perhaps not desired
    # df0 = discfunc.create("Adaptive",lagrangespace,name="test")
    # df0.clear();
    # df00 = discfunc.create("Adaptive",lagrangespace,name="testA")
    # df00.clear();
    # df0.addToVTKWriter(vtk,vtk.PointData)
    # df00.addToVTKWriter(vtk,vtk.PointData)

    expr_global = lambda x: [-(x[1] - 0.5)*math.sin(x[0]*12)]*dimRange

    print("HALLO C")
    gf = create.function("global", grid2d, "expr_global", 1, expr_global)
    print(gf)
    print(type(gf))
    print("***************************")
    print(dimRange,gf.__repr__(),len(expr_global([0,0])))
    print("***************************")
    print("HALLO D")
    df1 = lagrangespace.interpolate( gf, name="test" )
    print(dir(df1))
    print("HALLO E")
    df2 = lagrangespace.interpolate( [5]*dimRange, name="zero" )
    print(type(df2))
    print("HALLO F")
    df3 = lagrangespace.interpolate( df1, name="copy" )
    print("HALLO G")
    df4 = lagrangespace.interpolate( lambda x:
            [(x-[0.5,0.5]).infinity_norm,]*dimRange, name="radius")
    print("HALLO H")
    # lagrangespace=0
    grid2d.writeVTK("space_demo", pointdata=[gf,df1,df2,df3,df4])
    print("HALLO I")

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
    # print("P1Bubble(2)")
    # testSpace(grid2d, "P1Bubble",2,1)


print("ALUSimplex")
test("ALUSimplex")
print("YASPGRID A")
test("Yasp")
print("YASPGRID B")
test("Yasp")
print("END")
