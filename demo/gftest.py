from mpi4py import MPI

import math
import dune.fem

from dune.models.gridfunction import gridFunction

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc = dune.fem.create.space("Lagrange", grid, dimrange=2, polorder=2)

code = """value[0] = sin(x[0]);
    value[1] = cos(x[1]);
"""
func = gridFunction(grid, code).get()

def initial(x):
    r  = (x-[6,6]).two_norm
    return [ 0 if r>0.3 else 1, -0.5 ]
initial_gf = grid.globalGridFunction("initial", initial)
solution   = spc.interpolate(initial_gf, name="solution")

solution.interpolate( func )
