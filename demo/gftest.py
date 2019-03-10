from mpi4py import MPI
import math

import dune.fem
import ufl
import dune.ufl
import dune.create as create

factor = 10.

grid = create.grid("ALUConform", dune.grid.cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)

func1 = """
double s = sin(xGlobal[0]);
double c = cos(@const:fac*xGlobal[1]);
value[ 0 ] = s*s;
value[ 1 ] = s*c;
value[ 2 ] = c*c;
"""
func2 = """
double cx = cos(xGlobal[0])*@jac:test[0][0];
double cy = cos(@const:fac*xGlobal[1]);
double sx = sin(xGlobal[0]);
double sy = sin(@const:fac*xGlobal[1]);
value[ 0 ][ 0 ] = cx*sx*@gf:test[1];
value[ 0 ][ 1 ] = 0;
value[ 1 ][ 0 ] = cx*cy*@gf:test[0];
value[ 1 ][ 1 ] = -@const:fac*sx*sy;
value[ 2 ][ 0 ] = 0;
value[ 2 ][ 1 ] = -2.*@const:fac*cy*sy;
"""
code = { 'eval': func1, 'jac': func2 }

# coeffFunc = create.function("global", grid, "global_velocity", 1, lambda x: [1,2])
# func = create.function("cpp", grid, "code", 3, code, coefficients={"test": coeffFunc} )
# func.setConstant("fac", [factor])

uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 2, field="double")
x = ufl.SpatialCoordinate(ufl.triangle)
coeff = ufl.Coefficient(uflSpace)
const = dune.ufl.Constant(0)
c = ufl.cos(const*x[1])
s = ufl.sin(x[0])
expr = ufl.as_vector([ s*s*coeff[0], s*c, c*c ])
coeffFunc = create.function("global", grid, "global_velocity", 0, lambda x: [1,2])
funcUFL = create.function("ufl", grid, "ufl", 1, expr, coefficients={coeff: coeffFunc})
const.value = factor

func = funcUFL # needed since cpp function doesn't work

space = create.space("lagrange", grid, dimRange=3, order=2)
solution = space.interpolate(funcUFL, name="solution")

def expr_global(x):
    return [math.sin(x[0])**2, math.sin(x[0])*math.cos(factor*x[1]), math.cos(factor*x[1])**2]
control = create.function("global", grid, "expr_global", 3, expr_global)

def expr_local(en,x):
    y = en.geometry.toGlobal(x)
    return func.localFunction(en).evaluate(x) - control.localFunction(en).evaluate(x)
difference = create.function("local", grid, "difference", 3, expr_local )

grid.writeVTK("gftest", pointdata=[control,func,funcUFL,solution,difference])
