from mpi4py import MPI

import math
import dune.fem

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")

func1 = """value[0] = sin(xGlobal[0]);
    value[1] = cos(xGlobal[1]);
"""
func2 = """@dimrange=2
    for( int i = 0; i < dimDomain; ++i )
    {
        value[ 0 ][ i ] = 2 * xGlobal[ i ];
        for( int j = 0; j < dimDomain; ++j )
            value[ 0 ][ i ] *= (i == j ? 1.0 : xGlobal[ j ]*xGlobal[ j ]);
    }
"""
code = { 'eval': func1, 'jac': func2 }
func = grid.function("code", code=code)

solution = grid.interpolate(func, space="Lagrange", order=2, name="solution")

def expr_global(x):
    return [math.sin(x[0]), math.cos(x[1])]
control = grid.function("expr_global", globalExpr=expr_global)

def expr_local(en,x):
    y = en.geometry.position(x)
    return func.localFunction(en).evaluate(x) - control.localFunction(en).evaluate(x)
difference = grid.function( "difference", localExpr=expr_local )

grid.writeVTK("gftest", pointdata=[control,func,solution,difference])
