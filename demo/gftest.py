from mpi4py import MPI

import math
import dune.fem

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")

code = """value[0] = sin(xGlobal[0]);
    value[1] = cos(xGlobal[1]);
"""
func = grid.function("code", code=code)

solution   = grid.interpolate(func, space="Lagrange", order=2, name="solution")

def expr_global(x):
    return [math.sin(x[0]), math.cos(x[1])]
control = grid.function("expr_global", globalExpr=expr_global)

class LocalDiff:
    def __init__(self):
      self.dimR = 2
    def __call__(self,en,x):
        y = en.geometry.position(x)
        return func.localFunction(en).evaluate(x) - control.localFunction(en).evaluate(x)
difference = grid.function( "difference", localExpr=LocalDiff() )

grid.writeVTK("gftest", pointdata=[control,func,solution,difference])
