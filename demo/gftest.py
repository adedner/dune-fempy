from mpi4py import MPI

import math
import dune.fem
import ufl
import dune.ufl
import dune.fem.function as gf

factor = 10.

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")

func1 = """
double s = sin(xGlobal[0]);
double c = cos(factor*xGlobal[1]);
value[ 0 ] = s*s;
value[ 1 ] = s*c;
value[ 2 ] = c*c;
"""
func2 = """
double cx = cos(xGlobal[0]);
double cy = cos(factor*xGlobal[1]);
double sx = sin(xGlobal[0]);
double sy = sin(factor*xGlobal[1]);
value[ 0 ][ 0 ] = 2*cx*sx;
value[ 0 ][ 1 ] = 0;
value[ 1 ][ 0 ] = cx*cy;
value[ 1 ][ 1 ] = -factor*sx*sy;
value[ 2 ][ 0 ] = 0;
value[ 2 ][ 1 ] = -2.*factor*cy*sy;
"""
func1 = func1.replace("factor",str(factor)) # use some mechanism setConstant later...
func2 = func2.replace("factor",str(factor)) # use some mechanism setConstant later...
code = { 'eval': func1, 'jac': func2 }

dimR = 2
coeffFunc = grid.function("global_velocity", order=1, globalExpr=lambda x: [1,2])
func = grid.function("code", 3, code=code, coefficients={"test":coeffFunc}, constants={"const": dimR} )
exit(1)

uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 2, field="double")
x = ufl.SpatialCoordinate(ufl.triangle)
coeff = ufl.Coefficient(uflSpace)
const = ufl.Constant(ufl.triangle)
c = ufl.cos(const*x[1])
s = ufl.sin(x[0])
expr = ufl.as_vector([ s*s*coeff[0], s*c, c*c ])
coeffFunc = grid.function("global_velocity", order=0, globalExpr=lambda x: [1,2])
funcUFL = grid.function("ufl", order=1, ufl=expr, coefficients={coeff:coeffFunc})
funcUFL.setConstant(const,[factor])

solution = grid.interpolate(funcUFL, space="Lagrange", order=2, name="solution")

def expr_global(x):
    return [math.sin(x[0])**2, math.sin(x[0])*math.cos(factor*x[1]), math.cos(factor*x[1])**2]
control = grid.function("expr_global", order=3, globalExpr=expr_global)

def expr_local(en,x):
    y = en.geometry.position(x)
    return funcUFL.localFunction(en).evaluate(x) - control.localFunction(en).evaluate(x)
difference = grid.function( "difference", order=3, localExpr=expr_local )

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
