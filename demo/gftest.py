from mpi4py import MPI

import math
import dune.fem

from dune.models.gridfunction import importFunction

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

code = """
std::cout << "hello world" << std::endl;
"""
 # value[0] = sin(xglobal[0]);
 # value[1] = cosevaluate(xglobal[0]);
c = importFunction(grid, code).get()
print(dir(c))
c.evaluate()
