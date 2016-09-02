from mpi4py import MPI

import math
import dune.fem

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc = dune.fem.create.space("Lagrange", grid, dimrange=2, polorder=2)

code = """value[0] = sin(xGlobal[0]);
    value[1] = cos(xGlobal[1]);
"""
func = grid.localFunction(code)

solution   = spc.interpolate([0,0,0], name="solution")

def expr_global(x):
    return [math.sin(x[0]), math.cos(x[1])]
control = grid.globalGridFunction("expr_global", expr_global)

grid.writeVTK("control", pointdata=[control])

solution.interpolate( func )
grid.writeVTK("gftest", pointdata=[solution])
