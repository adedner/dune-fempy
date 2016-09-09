from mpi4py import MPI

import math
import dune.fem
import ufl
import dune.ufl
import dune.fem.function as gf

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")

func1 = """
double c = cos(xGlobal[1]);
double s = sin(xGlobal[0]);
value[ 0 ] = s*s;
value[ 1 ] = s*c;
value[ 2 ] = c*c;
"""
func2 = """double cx = cos(xGlobal[0]);
    double cy = cos(xGlobal[1]);
    double sx = sin(xGlobal[0]);
    double sy = sin(xGlobal[1]);
    value[ 0 ][ 0 ] = 2*cx*sx;
    value[ 0 ][ 1 ] = cx*cy;
    value[ 0 ][ 2 ] = 0;
    value[ 1 ][ 0 ] = 0;
    value[ 1 ][ 1 ] = -sx*sy;
    value[ 1 ][ 2 ] = -2*cy*sy;
    for( int j = 0; j < dimDomain; ++j )
        value[ 2 ][ j ] = 0;
"""
code = { 'eval': func1, 'jac': func2 }

dimR = 2
isConst = False
coef = { 'test': (dimR, isConst) }

coeffFunc = grid.function("global_velocity", globalExpr=lambda x: [1,2])
func = grid.function("code", coef, code=code)

# this is basically just a hack for now
class YClass( object ):
    pass
y= YClass()
setattr( y, 'number', 0 )
y.number

func.setCoefficient(y, coeffFunc)

uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 2, field="double")
x = ufl.SpatialCoordinate(ufl.triangle)
coeff = ufl.Coefficient(uflSpace)
c = ufl.cos(x[1])
s = ufl.sin(x[0])
expr = ufl.as_vector([ s*s*coeff[0], s*c, c*c ])
funcUFL = grid.function("ufl", ufl=expr)
funcUFL.setCoefficient(coeff, coeffFunc)

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
