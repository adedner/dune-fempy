
# coding: utf-8

# # Grid functions (and some plotting) [(Notebook)][1]
#
# [1]: _downloads/gridfunctions.ipynb
# *Grid function* are perhaps the most fundamental concept of (after the grid itself). These are *localizable* function, i.e., functions that can be evaluated given an element `e` of the grid and a coordinate in the reference element of `e`. To this end they provide a method `lf = localfunction(e)` returning an object with an `evaluate(x),jacobian(x),hessian(x)` method. The coordinate `x` should be in the reference element of `e` but note that the derivates returned are with respect to the coordinate system of `e` itself. All numerical algorithm provided in the dune modules do not use globally defined functions but rely on *grid function*.
#
# In addition we also show some simple inline plotting functionality.

# In[ ]:


try:
    get_ipython().magic('matplotlib inline # can also use notebook or inline # can also use notebook or nbagg')
except:
    pass
import math
from dune.grid import gridFunction
import dune.fem
import ufl
import dune.ufl
import dune.create as create
from dune.fem.plotting import plotPointData as plot
from dune.fem.plotting import plotComponents as plotComponents
from dune.fem.space import lagrange as lagrangeSpace
from dune.fem.function import globalFunction, localFunction, uflFunction, cppFunction

grid = create.grid("ALUConform", dune.grid.cartesianDomain([0,0],[1,1],[4,4]), dimgrid=2)


# There are a number of different ways *grid function* can be constructed. We start with the simplest approach using callback to python function:

# In[ ]:


def expr_global(x):
    return [math.sin(x[0])**2, math.sin(x[0])*math.cos(10*x[1]), math.cos(10*x[1])**2]
gf_gl = globalFunction(grid, "expr_global", 3, expr_global)
print(expr_global([1,1]),gf_gl(grid.elements.__next__(),[0.5,0.5]))
# plot all components + the grid
plotComponents(gf_gl, gridLines="black")
# plot the absolute value of the compnents without the grid
plot(gf_gl,gridLines="")

def expr_local(en,x):
    y = en.geometry.toGlobal(x)
    return [gf_gl.localFunction(en).evaluate(x)[1] - y.two_norm,10]
gf_loc = localFunction( grid, "expr_local", 3, expr_local)
# here we show the solution (first component) overlayed with the grid (default)
plot(gf_loc[0])


# Note how we can access the one component of a grid function by writing `gf_loc[0]`.
#
# This is quite a coarse grid. By changing the `level` argument in the plotting functions one can get a better representation of the functions without changing the grid. Since the data was never interpolate onto a grid the full information of the original function is still available:

# In[ ]:


plot(gf_loc[0],level=3)
# of coarse we could refine the grid and plot with the original level
grid.hierarchicalGrid.globalRefine(4)
plot(gf_loc[0],level=3)


# A few remarks: note that all grid functions return a vector - even scalar functions as the ones in the example. The dimension of the range of these functions is automatically deduced from its return value.
#
# The constructor arguments for grid functions are always similar: in addition to the grid, a name is used (mostly for visualization purposes), a polynomial order is used (for detecting quadrature orders in numerical schemes); the last argument the the expression defining the actual function.
#
# Finaly, since evaluating these functions requires a callback from C++ to python and no *just in time* compilation is used, this approach can be quite slow compared to some of the others described in the following.

# In[ ]:


uflSpace = dune.ufl.Space(grid, 2, field="double")
x = ufl.SpatialCoordinate(ufl.triangle)
const = dune.ufl.Constant(0)
c = ufl.cos(const*x[1])
s = ufl.sin(x[0])
expr = ufl.as_vector([ s*c ])
funcUFL = uflFunction(grid, "ufl", 1, expr)
const.value = 10
# in the case of a scalar function the actual values are plotted
# and not the absolute values are for dimRange > 1
plot(funcUFL)

scalarUFLSpace = dune.ufl.Space(grid,1)
coeff = ufl.Coefficient(scalarUFLSpace)
coeffFunc = globalFunction(grid, "global_velocity", 0, lambda x: [2])
expr = ufl.as_vector([ s*s*coeff[0] ])
funcUFL1 = uflFunction(grid, "ufl1", 4, expr,
#                           coefficients={coeff: funcUFL})
                           coefficients={coeff: coeffFunc})
plot(funcUFL1)
# for e in grid.elements:
#     print( funcUFL1.localFunction(e) )


# Finally the grid function can be described proving strings with valid C++ code

# In[ ]:


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

# Not working at the moment
# coeffFunc = globalFunction(grid, "global_velocity", 1, lambda x: [1,2])
# func = cppFunction(grid, "code", 3, code, coefficients={"test": coeffFunc} )
# func.setConstant("fac", 10)
# show all components but not the grid
# plotComponents(func,gridLines="")


# Most of the methods available on the grid functions have been used above - mainly the `localfunction` method. There are a few more of interst are:

# In[ ]:


# print("dimension of range: ",func.dimRange)
# print("underlying grid instance: ",func.grid)
# print("name: ",func.name)


# The following two methods can be used together with the `grid.triangulation( level )` method to obtain a simple numpy vector representation of the solution:

# In[ ]:


# print("return array with values at cell centers: ", func.cellData(0))
# print("return array with values at nodes: ", func.pointData(0))


# The following method computes the integral of the function over the domain - a possible usage is to compute the approximation error of a function:

# In[ ]:


# print("value of integral over whole domain: ", func.integrate())


# The above results in grid functions which are still identical to the global functions that they represent, i.e., no discretization is involved. The next examples use an underlying discrete function space and interpolates either global functions or grid functions into this space:

# In[ ]:


# first a constant function
vector_space = lagrangeSpace(grid, dimRange=grid.dimension, order=2)
constant = vector_space.interpolate([1,1], name="constant")
uh = vector_space.interpolate(lambda x: [math.sin(x[0]*x[1]*math.pi),math.cos(x[0]*x[1]*math.pi)], name="xy")
plot(uh,vectors=[0,1],gridLines="")
scalar_space = lagrangeSpace(grid, dimRange=1, order=2)
# we can also use the grid function decorators to define the a grid function
vorticity = gridFunction(scalar_space.grid)(
            lambda en,x: [uh.localFunction(en).jacobian(x)[1][0]-uh.localFunction(en).jacobian(x)[0][0]])
vorticity_h = scalar_space.interpolate(vorticity, name="fh")
plot(vorticity_h)


# Note the difference between the grid function `vorticity` and `vorticity_h` - the second being an interpolation of the former into the discrete function space.
#
# While discrete functions are grid functions and have the same attributes, they have a few additional methods

# In[ ]:


print("number of dofs: ", vorticity_h.size)
print("the underlying space: ", vorticity_h.space)
# we can construct a new discrete function of the same type and the same values
copy = vorticity_h.copy()
# set all value to zero
copy.clear()
# assign all dofs to be identical to another discrete function over the same space
copy.assign(vorticity_h)
# store the interpolation of some (grid) function
copy.interpolate(vorticity_h)  # the oritignal discrete vorticity
copy.interpolate(vorticity)    # does the same by interpolating the original grid function
copy.interpolate([1])          # copy(x) = 1
copy.interpolate(lambda x: [x[0]*x[1]] )     # interpolation of a global function
copy.interpolate(lambda en,x: [en.geometry.center.two_norm] ) # interpolation of a elementwise function
plot(copy)


# Note that in the last case a piecewise constant function is interpolated onto the continuous first order lagrange space.
#
# Underlying any discrete function is a vector containing the degrees of freedom. This can also retrieved using the `dofVector` property on the discrete function. In general this will not be of much help since there are no access methods provided at the moment.
#
# Accessing the dof vector becomes useful when the `fem` or `eigen` backend is used for storing the dofs. In this case the
# python buffer protocol is used to convert the dof vector into an `numpy` array without requiring any copy.
# This requires having the [Eigen package][1].
#
# [1]: http://eigen.tuxfamily.org/index.php?title=Main_Page

# In[ ]:


import numpy as np
spc = lagrangeSpace(grid, dimRange=1, order=1, storage='fem')
uh = spc.interpolate(vorticity_h)
plot(uh)
uh_dofs = uh.as_numpy # equivalent to np.array( uh.dofVector, copy=False )
uh_dofs *= -2
plot(uh)


# We can also use an existing `numpy` array as basis of a discrete function (over space with `fem` or `eigen` backend) - this makes it easy to use for example `scipy` solvers directly:

# In[ ]:


# a simple function to check the address of the memory block of a numpy array
def id(x):
    return x.__array_interface__['data'][0]

# set up a numpy vector with random entries
dofs = np.random.random(spc.size)
print("id of numpy array: ", id(dofs))
# reinterpret the numpy vector as a discrete function
xh = spc.function("name", dofVector=dofs)
# get the dof vector as numpy array
xh_dofs = xh.as_numpy
print("id of dof vector: ", id(xh_dofs))

# plot the random values
plot(xh)
# set the dof ventries unsing the discrete vorticity
xh.interpolate(vorticity_h)
plot(xh)
# reset all entries in the dof vector to random values
dofs[:] = np.random.random(spc.size)
plot(xh) # constructor of numpy array called here, why?
# now do the same using the extracted numpy array
xh_dofs[:] = np.random.random(spc.size)
plot(xh)
