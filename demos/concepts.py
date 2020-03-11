# <markdowncell>
# # Setting up the Mesh
# After some general import statements we start our tutorial by setting up
# a simple Cartesian grid over the domain $[0,1]^2$ subdivided into four
# intervals in each dimension. We can then show the grid using the `plot`
# method - note that direct plotting with MatPlotLib is only supported in
# 2D, other options for visualization, e.g., using ParaView will be
# discussed later.
#
# Using the `dune.grid.structuredGrid` function, is the simplest way to construct
# a grid, more complicated
# (unstructured) grids can be constructed using a dictionary containing the
# vertices and the element connectivity or by for example providing a
# `gmsh` file. This will be discussed in later section, see for example
# here for [a 3d gmsh file](othergrids_nb.ipynb#A-3D-Example-Using-a-GMesh-File)
# or here [for a simple example using a dictonary](dune-corepy_nb.ipynb#An-unstructured-grid)

# <codecell>

import time, numpy, math, sys, io, matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot

from dune.grid import structuredGrid as leafGridView
gridView = leafGridView([0, 0], [1, 1], [4, 4])
gridView.plot()

# <markdowncell>
# # Grid Functions
# These are function that are defined over a given grid and are evaluated
# by using an element of the grid and local coordinated within that
# element:
# ```
#    value = gridFunction(element,localCoordinate)
# ```
# Alternatively one can obtain a `LocalFunction` from a grid function which
# can be bound to an element and then evaluate via local coordinate:
# ```
#    localFunction = gridFunction.localFunction()
#    for e in grid.elements:
#        localFunction.bind(e)
#        value = localFunction(x)
# ```
# There are multiple ways to construct grid functions. The easiest way it
# to use UFL expression. Many methods expecting grid functions are argument
# can also directly handle UFL expression.
# We can for example integrate a UFL expression over the grid:
# <codecell>

from ufl import SpatialCoordinate, triangle
x = SpatialCoordinate(triangle)

exact = 1/2*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1

from dune.fem.function import integrate
print( integrate(gridView, exact, order=5) )

# <markdowncell>
# and plot them using matplotlib or write a vtk file for postprocessing
# <codecell>

from dune.fem.plotting import plotPointData as plot
plot(exact, grid=gridView)
gridView.writeVTK('exact', pointdata={'exact': exact})

# <markdowncell>
# In some cases it can be necessary to convert a UFL expression into a grid
# function explicitly - for example to be able to evaluate it over each
# element in the grid:
# <codecell>

from dune.fem.function import uflFunction
exact_gf = uflFunction(gridView, name="ufl", order=1, ufl=exact)
mass = 0
for element in gridView.elements:
    mass += exact_gf(element,[0.5,0.5]) * element.geometry.volume
print(mass)

# <markdowncell>
# Another way to obtain a grid function is to use the `gridFunction`
# decorator. This can be obtained from `dune.grid` but then without UFL
# support. Using the decorator from `dune.fem.function` the resulting grid
# function can be used seamlessly within UFL expressions:
# <codecell>

from dune.fem.function import gridFunction
@gridFunction(gridView,name="callback",order=1)
def exactLocal(element,xLocal):
    x = element.geometry.toGlobal(xLocal)
    return 1/2.*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1

# <markdowncell>
# we can use the same approach but with a function using global
# coordinates but can then be used like any other grid function:
# <codecell>

@gridFunction(gridView,name="callback",order=1)
def exactGlobal(x):
    return 1/2.*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1

lf = exactGlobal.localFunction()
mass = 0
for element in gridView.elements:
    lf.bind(element)
    mass += lf([0.5,0.5]) * element.geometry.volume
print(mass)

print( integrate(gridView, [exact,exactLocal,exactGlobal], order=5) )

# <markdowncell>
# As pointed out the `dune.fem` grid function can be used like any other
# UFL coefficient to form UFL expressions:
# <codecell>

print( integrate(gridView, abs(exact-exactLocal), order=5) )
gf = uflFunction(gridView, name="ufl", order=1, ufl=exact+exactLocal*exactGlobal)
fig = pyplot.figure(figsize=(20,10))
gf.plot(figure=(fig,121))
exactGlobal.plot(figure=(fig,122))

# <markdowncell>
# Converting UFL expressions to grid functions leads to JIT code generation
# and is therefore efficient when used in other C++ algorithm (like
# `integrate`). On the other hand using the `gridFunction` decorator leads
# to a callback into Python for each evaluation and is therefore much less
# efficient. An alternative approach is based on writing small C++ snippets
# implementing the grid function:
# <codecell>

from dune.fem.function import cppFunction
code="""
#include <cmath>
#include <dune/common/fvector.hh>
template <class GridView>
auto aTimesExact(double a) {
  // the return value needs to be a lambda which either returns a `double`
  // or a `Dune::FieldVector<double,R>` with `R>=1`
  return [a](const auto& en,const auto& xLocal) -> auto {
    auto x = en.geometry().global(xLocal);
    return a*(1./2.*(std::pow(x[0],2)+std::pow(x[1],2)) - 1./3.*(std::pow(x[0],3) - std::pow(x[1],3)) + 1.);
  };
}
"""
exactCpp = cppFunction(gridView, name="exactCpp", order=2,
                       fctName="aTimesExact",includes=io.StringIO(code))(2.)
print( integrate(gridView, abs(2*exact-exactCpp), order=5) )
exactCpp.plot()

# <markdowncell>
# As the above example shows it is easy to pass in parameters to the C++
# implementation - here a double `2`. Note that it is not possible to
# obtain a reference to this double so to make sure a change to the
# constant on the Python side carries over to the C++ side an option is to
# use a `dune.common.FieldVector` instead.
# These parameters can also be more complex e.g. other grid # function -
# note that a UFL expression can not be directly passed in - it first needs
# to be converted into a grid function using `uflFunction`.
# <codecell>

code2="""
#include <cmath>
#include <dune/common/fvector.hh>
template <class GridView, class GF>
auto aTimesExact(const GF &gf,Dune::FieldVector<double,1> &a) {
  return [lgf=localFunction(gf),&a] (const auto& en,const auto& xLocal) mutable -> auto {
    lgf.bind(en); // lambda must be mutable so that the non const function can be called
    return a[0]*lgf(xLocal);
  };
}
"""
from dune.common import FieldVector
a = FieldVector([2])
exactCpp2 = cppFunction(gridView, name="exactCpp", order=2,
                        fctName="aTimesExact",includes=io.StringIO(code2))(exact_gf,a)
print( integrate(gridView, abs(exactCpp-exactCpp2), order=5) )
a[0] = 4
print( integrate(gridView, abs(2*exactCpp-exactCpp2), order=5) )

# <markdowncell>
# The above is just one of a few ways C++ code snippets can be used
# together with Python code to improve efficiency or extend the existing
# binding to Dune. In the above example there is no advantage of using the
# C++ code over the code generated based on the UFL expression. For more
# complicated functions e.g. with many if statements or based on more
# information from the given element like it's neighbors the expressibility
# of UFL might not be sufficient or lead to hard to read code. In these
# cases directly providing C++ code (or Python code) can be a reasonable
# alternative.
# <codecell>

# <markdowncell>
# # Discrete Function Spaces
# Note that the grid functions set up so far did not involve any
# discretization, they are exactly evaluated at the given point.
# A special type of grid functions are discrete functions living in a
# discrete (finite dimensional) space.
# <codecell>

from dune.fem.space import lagrange as solutionSpace
space = solutionSpace(gridView, order=2)

# <markdowncell>
# The easiest way to construct a discrete function is to use the interpolate
# method on the discrete function space.
# <codecell>

u_h = space.interpolate(exact, name='u_h')

# <markdowncell>
# On an existing discrete function the `interpolate` method can be used to
# reinitialize it
# <codecell>

u_h.interpolate( cppFunction(gridView, name="exactCpp", order=2,
                 fctName="aTimesExact",includes=io.StringIO(code))(2.) )
u_h.interpolate( lambda x: 1/2.*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1 )

# <markdowncell>
# Note that in the last example the Python lambda is used as a callback
# automatically using the same concept used in the `gridFunction`
# decorator. As pointed out above there are some methods where these
# conversions are implicit and no explicit generation of a grid function
# has be carried out.
#
# If a discrete function is already available it is possible to call `copy`
# to obtain further discrete functions:
# <codecell>

u_h_n = u_h.copy(name="previous")

# <markdowncell>
# Finally, `clear` can be called on a discrete function which sets all
# coefficient to zero and `assign` can be used to copy all coefficients
# between two discrete function over the same space:
# <codecell>

u_h_n.clear()
u_h_n.assign( u_h )

# <markdowncell>
# All the things we did above with grid functions can be done with discrete
# functions, e.g., evaluate locally
# <codecell>

localUh = u_h.localFunction()
mass = 0
for element in gridView.elements:
    localUh.bind(element) # using u_h(element,[0.5,0.5]) also works
    mass += localUh([0.5,0.5]) * element.geometry.volume
print(mass)

# <markdowncell>
# or plot using matplotlib and write a vtk file for postprocessing
# <codecell>

u_h.plot(gridLines="white")
gridView.writeVTK('uh', pointdata=[u_h])

# <markdowncell>
# and as before a discrete function can be used as coefficient in a UFL
# expression
# <codecell>

print( integrate(gridView, abs(exact-u_h), order=5) )

# <markdowncell>
# Note: the discrete function `u_h` already has a `name` attribute given in
# the `interpolate` call. This is used by default in the vtk file. An
# alternative name can be given by using a dictionary as shown previously.
# <codecell>

# <markdowncell>
# # Models and Schemes
# We consider a scalar boundary value problem
# \begin{align*}
# -\triangle u &= f & \text{in}\;\Omega:=(0,1)^2 \\
# \nabla u\cdot n &= g_N & \text{on}\;\Gamma_N \\
# u &= g_D & \text{on}\;\Gamma_D
# \end{align*}
# and $f=f(x)$ is some forcing term.
# For the boundary conditions we set $\Gamma_D={0}\times[0,1]$ and take
# $\Gamma_N$ to be the remaining boundary of $\Omega$.
#
# We will solve this problem in variational form
# \begin{align*}
# \int \nabla u \cdot \nabla \varphi \
# - \int_{\Omega} f(x) \varphi\ dx
# - \int_{\Gamma_N} g_N(x) v\ ds
# = 0.
# \end{align*}
# We choose $f,g_N,g_D$ so that the exact solution
# is given by
# \begin{align*}
# u(x) = \left(\frac{1}{2}(x^2 + y^2) - \frac{1}{3}(x^3 - y^3)\right) + 1
# \end{align*}
# <codecell>

from ufl import TestFunction, TrialFunction
from dune.ufl import DirichletBC
u = TrialFunction(space)
v = TestFunction(space)

from ufl import dx, grad, div, grad, dot, inner, sqrt, conditional, FacetNormal, ds
a = dot(grad(u), grad(v)) * dx

f   = -div( grad(exact) )
g_N = grad(exact)
n   = FacetNormal(space)
b   = f*v*dx + dot(g_N,n)*conditional(x[0]>=1e-8,1,0)*v*ds
dbc = DirichletBC(space,exact,x[0]<=1e-8)

# <markdowncell>
# With the model described as a ufl form, we can construct a scheme class
# that provides the solve method which we can use to compute the solution:
# <codecell>

from dune.fem.scheme import galerkin as solutionScheme
scheme = solutionScheme([a == b, dbc], solver='cg')
scheme.solve(target = u_h)

# <markdowncell>
# We can compute the error between the exact and the discrete solution by
# using the `integrate` function described above:
# <codecell>

h1error = dot(grad(u_h - exact), grad(u_h - exact))
error = sqrt(integrate(gridView, h1error, order=5))
print("Number of elements:",gridView.size(0),
      "number of dofs:",space.size,"H^1 error:", error)

# <markdowncell>
# To verify that the discrete scheme is converging to the exact solution
# we can compute the experimental order of convergence (EOC):
# \begin{align*}
# {\rm eoc} = \frac{\log{e_h/e_H}}{\log{h/H}}
# \end{align*}
# where $h,H$ refer to the spacing of two grids and $e_h,e_H$ are measures
# for the accuracy of the discrete solution, i.e., $e_h=\|u-u_h\|,e_H=\|u-u_h$
# using a suitable norm (here $h_1$ as used before):
# <codecell>

from math import log
loops = 2
for eocLoop in range(loops):
    error_old = error
    gridView.hierarchicalGrid.globalRefine(1)
    u_h.interpolate(0)
    scheme.solve(target = u_h)
    error = sqrt(integrate(gridView, h1error, order=5))
    eoc = round(log(error/error_old)/log(0.5),2)
    print("EOC:",eoc,
          "Number of elements:", gridView.size(0),
          "number of dofs:", space.size,"H^1 error:", error)

# <markdowncell>
# We have already seen how grid function can be used within UFL expressions
# and forms - in the next section we will give another example for this in the
# context of a time dependent problems. In addition we can also use the
# `Constant` class to add constants to UFL expressions/forms which then can
# be changed without requiring any recompilation of the model. An example
# would again be in a time dependent problem a time varying coefficient.
# Being able to change the value of the time in the model without
# recompilation is crucial for an efficient code. We will demonstrate this
# here by adding a mass term with a constant mass $m$. We will not change
# the right hand side so the exact solution we used so far will only be
# valid for $m=0$.
# <codecell>

from dune.ufl import Constant

m = Constant(1,name="mass")   # start with m=1
a = dot(grad(u), grad(v)) * dx + m*u*v * dx
scheme = solutionScheme([a == b, dbc], solver='cg')

scheme.solve(target = u_h)

h1error = dot(grad(u_h - exact), grad(u_h - exact))
error = sqrt(integrate(gridView, h1error, order=5))
print("Number of elements:",gridView.size(0),
      "number of dofs:",space.size,"H^1 error:", error)

# <markdowncell>
# We can print the value of a `Constant` with name `foo` via
# `scheme.model.foo` and change it's value using the same attribute:
# <codecell>

print(scheme.model.mass)
scheme.model.mass = 0      # go back to original problem
print(scheme.model.mass)

scheme.solve(target = u_h)

h1error = dot(grad(u_h - exact), grad(u_h - exact))
error = sqrt(integrate(gridView, h1error, order=5))
print("Number of elements:",gridView.size(0),
      "number of dofs:",space.size,"H^1 error:", error)

# <markdowncell>
# # Listing Available Dune Components
# The available realization of a given interface, i.e., the available
# grid implementations, depends on the modules found during configuration.
# Getting access to all available components is straightforward:
# <codecell>

from dune.utility import components
# to get a list of all available components:
components()
# to get for example all available grid implementations:
components("grid")
