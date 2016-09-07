from __future__ import print_function
import math
from mpi4py import MPI

import dune.fem as fem

def testSpace(gridtype):
    grid2d = fem.leafGrid("../data/unitcube-2d.dgf", gridtype, dimgrid=2)

    def expr_global(x):
        return [-(x[1] - 0.5)*math.sin(x[0]*12),x[0]*x[1]]

    gf  = grid2d.function("expr_global", order=1, exprGlobal=expr_global)
    df1 = grid2d.interpolate(gf, space="Lagrange",name="interpolate")
    df2 = grid2d.interpolate([5,3], space="Lagrange") # , storage="Numpy" ) # , name="53" )
    df3 = grid2d.interpolate(df1, space="Lagrange", name="copy", storage="Istl" )
    df4 = grid2d.interpolate(lambda x: [expr_global(x)[0]], space="Lagrange", name="test")
    df4 = grid2d.interpolate(gf, space="Lagrange", name="test2", polorder=2)

    grid2d.writeVTK("space_demo", pointdata=[gf,df1,df2,df3,df4])

print("ALUGRID")
testSpace("ALUSimplexGrid")
print("YASPGRID A")
testSpace("YaspGrid")
print("YASPGRID B")
testSpace("YaspGrid")
print("END")
