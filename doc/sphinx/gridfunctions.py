# coding: utf-8

# ## Grid functions (and some plotting)
#
# *Grid function* are perhaps the most fundamental concept of (after the grid itself). These are *localizable* function, i.e., functions that can be evaluated given an element `e` of the grid and a coordinate in the reference element of `e`. To this end they provide a method `lf = localfunction(e)` returning an object with an `evaluate(x),jacobian(x),hessian(x)` method. The coordinate `x` should be in the reference element of `e` but note that the derivates returned are with respect to the coordinate system of `e` itself. All numerical algorithm provided in the dune modules do not use globally defined functions but rely on *grid function*.
#
# In addition we also show some simple inline plotting functionality.

# In[1]:

import math
import dune.fem
import ufl
import dune.ufl
import dune.create as create
from dune.fem.ipython import plotPointData as plot
from dune.fem.ipython import plotComponents as plotc

grid = create.grid("ALUConform", dune.grid.cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)


# There are a number of different ways *grid function* can be constructed. We start with the simplest approach using callback to python function:

# In[2]:

def expr_global(x):
    return [math.sin(x[0])**2, math.sin(x[0])*math.cos(10*x[1]), math.cos(10*x[1])**2]
gf_gl = create.function("global", grid, "expr_global", 3, expr_global)
# plot all components + the grid
plotc(grid, gf_gl)
# plot the absolute value of the compnents without the grid
plot(grid,gf_gl,component=-1,showGrid=False)

def expr_local(en,x):
    y = en.geometry.position(x)
    return [gf_gl.localFunction(en).evaluate(x)[1] - y.two_norm,10]
gf_loc = create.function("local", grid, "expr_local", 3, expr_local)
# here we show the solution (first component is default) overlayed with the grid
plot(grid,gf_loc)


# A few remarks: note that all grid functions return a vector - even scalar functions as the ones in the example. The dimension of the range of these functions is automatically deduced from its return value.
#
# The constructor arguments for grid functions are always similar: in addition to the grid, a name is used (mostly for visualization purposes), a polynomial order is used (for detecting quadrature orders in numerical schemes); the last argument the the expression defining the actual function.
#
# Finaly, since evaluating these functions requires a callback from C++ to python and no *just in time* compilation is used, this approach can be quite slow compared to some of the others described in the following.

# In[3]:

uflSpace = dune.ufl.Space(grid, 2, field="double")
x = ufl.SpatialCoordinate(ufl.triangle)
coeff = ufl.Coefficient(uflSpace)
const = ufl.Constant(ufl.triangle)
c = ufl.cos(const*x[1])
s = ufl.sin(x[0])
expr = ufl.as_vector([ s*c])
funcUFL = create.function("ufl", grid, "ufl", 1, expr)
funcUFL.setConstant(const, [10])
plot(grid, funcUFL)

coeffFunc = create.function("global", grid, "global_velocity", 0, lambda x: [1,2])
expr = ufl.as_vector([ s*s*coeff[0] ])
funcUFL1 = create.function("ufl", grid, "ufl1", 1, expr,
                           #coefficients={coeff: funcUFL})
                           coefficients={coeff: coeffFunc})
plot(grid, funcUFL1)


# The above results in grid functions which are still identical to the global functions that they represent, i.e., no discretization is involved. The next examples use an underlying discrete function space and interpolates either global functions or grid functions into this space:

# In[4]:

# first a constant function
vector_space = create.space("Lagrange", grid, dimrange=grid.dimension, order=2)
constant = vector_space.interpolate([1,1], name="constant")
uh = vector_space.interpolate(lambda x: [math.sin(x[0]*x[1]*math.pi),math.cos(x[0]*x[1]*math.pi)], name="xy")
plot(grid,uh,component=-2,showGrid=False)
scalar_space = create.space("Lagrange", grid, dimrange=1, order=2)
vorticity = create.function("local", grid, "global_velocity", 3,
                       lambda en,x: [uh.localFunction(en).jacobian(x)[1][0]-uh.localFunction(en).jacobian(x)[0][0]])
vorticity_h = scalar_space.interpolate(vorticity, name="fh")
plot(grid,vorticity_h)


# Finally the grid function can be described proving strings with valid C++ code

# In[5]:

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

coeffFunc = create.function("global", grid, "global_velocity", 1, lambda x: [1,2])
func = create.function("cpp", grid, "code", 3, code, coefficients={"test": coeffFunc} )
func.setConstant("fac", [10])
plot(grid,func)


# In[ ]:
