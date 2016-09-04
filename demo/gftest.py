from mpi4py import MPI

import math
import dune.fem
import ufl

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")

code = """
    double c = cos(xGlobal[1]);
    double s = sin(xGlobal[0]);
    value[ 0 ] = s*s;
    value[ 1 ] = s*c;
    value[ 2 ] = c*c;
"""
func = grid.function("code", code=code)

x = ufl.SpatialCoordinate(ufl.triangle)
c = ufl.cos(x[1])
s = ufl.sin(x[0])
expr = ufl.as_vector([ s*s, s*c, c*c ])
funcUFL = grid.function("ufl", ufl=expr)

solution = grid.interpolate(func, space="Lagrange", order=2, name="solution")

def expr_global(x):
    return [math.sin(x[0])**2, math.sin(x[0])*math.cos(x[1]), math.cos(x[1])**2]
control = grid.function("expr_global", globalExpr=expr_global)

def expr_local(en,x):
    y = en.geometry.position(x)
    return func.localFunction(en).evaluate(x) - control.localFunction(en).evaluate(x)
difference = grid.function( "difference", localExpr=expr_local )

# method 1
grid.writeVTK("gftest", pointdata=[control,func,funcUFL,solution,difference])

# method 2
vtk = grid.vtkWriter()
grid = None   # is not needed anymore on python side
control.addToVTKWriter(vtk, dune.common.DataType.CellData)
func.addToVTKWriter(vtk, dune.common.DataType.CellData)
funcUFL.addToVTKWriter(vtk, dune.common.DataType.CellData)
solution.addToVTKWriter(vtk, dune.common.DataType.CellData)
difference.addToVTKWriter(vtk, dune.common.DataType.CellData)
vtk.write( "gftest" )
