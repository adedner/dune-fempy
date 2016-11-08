# coding: utf-8

# # Solving an Elliptic PDE
# ## Getting started
#
# This demo introduces basic usage of dune-fempy, using the Poisson equation as an example. Namely,
#
# \begin{align*}
#   - \Delta u + u &= f && \text{in } \Omega \\
#   \nabla u \cdot \textbf{n} &= 0 && \text{on } \Gamma
# \end{align*}
#
#
# If you have compiled DUNE against MPI, we strongly advise you to first initialize MPI from Python.
# At least OpenMPI is known to fail, if initialized only in the dune-fempy library.

# In[1]:

import dune.fem
dune.fem.parameter.append("../data/parameter")


# First, we create our computational grid. Our domain will be the unit square divided into 8x8 quadrilaterals. To actually create the grid, we choose an implementation of the DUNE grid interface: a 2-dimensional ALUGrid with simplices and conforming bisection refinement.

# In[2]:

import dune.create as create
grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8]), dimgrid=2)


# We set up the base variables u, v and x in UFL.

# In[3]:

from dune.ufl import Space
from ufl import TestFunction, TrialFunction, SpatialCoordinate
uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())


# Next we define the equation for the weak form, given by
#
# \begin{equation}
# \int_{\Omega} uv + \nabla u\cdot\nabla v \ dx =  \int_{\Omega} f v \ dx.
# \end{equation}
# We take $f = 9\pi^2\cos(2\pi x_0)\cos(2\pi x_1)$.
#
# Note that then the exact solution is then
# $u = \cos(2\pi x_0)\cos(2\pi x_1)$.

# In[4]:

from math import pi
from ufl import cos, as_vector, dx, grad, inner
f = 9*pi*pi*cos(2*pi*x[0])*cos(2*pi*x[1])
exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )
equation = (inner(grad(u), grad(v)) + inner(u,v)) * dx == f * v[0] * dx
equation


# We create the space and the model.

# In[5]:

spc = create.space("Lagrange", grid, dimrange=1, order=1)
model = create.model("elliptic", grid, equation)


# We create the scheme and set parameters for the solver.

# In[6]:

scheme = create.scheme("h1", spc, model,       parameters=       {"fem.solver.newton.linabstol": 1e-10,
        "fem.solver.newton.linreduction": 1e-10,
        "fem.solver.newton.verbose": 0,
        "fem.solver.newton.linear.verbose": 0})


# We create a grid function for our exact solution.

# In[7]:

exact_gf = create.function("ufl", grid, "exact", 5, exact)


# We set up a function for plotting the data using matplotlib.

# In[8]:

try:
    import matplotlib
    from matplotlib import pyplot
    from numpy import amin, amax, linspace
    from IPython.core.display import display

    def plot(grid, solution):
        triangulation = grid.triangulation(4)
        data = solution.pointData(4)

        levels = linspace(amin(data[:,0]), amax(data[:,0]), 256)

        fig = pyplot.figure()
        fig.gca().set_aspect('equal')
        pyplot.triplot(grid.triangulation(), antialiased=True, linewidth=0.2, color='black')
        pyplot.tricontourf(triangulation, data[:,0], cmap=pyplot.cm.rainbow, levels=levels)
        display(pyplot.gcf())
except ImportError as e:
    print(e)
    def plot(grid, solution):
        pass


# Now we solve the system. We assign the solution to `uh`, and define a function to calculate the $L^2$ error, i.e. $|u_h - u|_{L^2}$. We output the data to a vtk file with name `laplace`, and plot it using `plot`. Finally we refine the grid twice and repeat the process.

# In[9]:

from math import sqrt
for i in range(2):
    print("solve on level", i, "number of dofs=", grid.size(2))
    uh,_ = scheme.solve()
    def l2error(en,x):
        val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
        return [ val[0]*val[0] ];
    l2error_gf = create.function("local", grid, "error", 5, l2error)
    error = sqrt(l2error_gf.integrate()[0])

    print("size:", grid.size(0), "L2-error:", error)
    grid.writeVTK("laplace", pointdata=[uh, l2error_gf])

    plot(grid, uh)

    if i < 1:
        grid.hierarchicalGrid.globalRefine(2)


# Congratulations! You have successfully solved and visualized your first PDE using dune-fempy.
#
# ## Different storage backends - using numpy/scipy
#
# For the first example we used solvers available in dune-fem - simple Krylov solvers with only diagonal preconditioning. Changing the `storage` argument in the construction of the space makes it possible to use more sophisticated solvers (either better preconditioners or direct solvers). For example
# ~~~
# spc = create.space("Lagrange", grid, dimrange=1, order=1, storage="istl")
# ~~~
# in the above code will switch to the solvers from `dune-istl`, other options are for example `eigen` or `petsc`.
#
# It is also possible to store the degrees of freedom in such a way that they can be treated as `numpy` vectors and an assembled system matrix can be stored in a `sympy` sparse matrix.
#
# __Note__: at the moment we can only export `Eigen` matrices to python so to get the following example to run, the `Eigen` package must be available and `dune-py` must have been configured with `Eigen`.
#
# Since we will be implementing a Newton solver first, let's study a truelly non linear problem - a version of the p-Laplace problem:
# \begin{gather}
#   - \frac{d}{2}\nabla\cdot |\nabla u|^{p-2}\nabla u + u = f
# \end{gather}

# In[10]:

import numpy as np
import scipy.sparse.linalg
import scipy.optimize
from ufl import ds

try:
    spc = create.space("Lagrange", grid, dimrange=1, order=1, storage='eigen')
except Exception as ex:
    template = "An exception of type {0} occured. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print(message)
    exit(1)
except:
    print('something...')
    exit(1)

d = 0.001
p = 1.7

rhs = (x[0] + x[1]) * v[0]
a = (pow(d + inner(grad(u), grad(u)), (p-2)/2)*inner(grad(u), grad(v)) + inner(u, v)) * dx + 10*inner(u, v) * ds
b = rhs * dx + 10*rhs * ds
model = create.model("elliptic", grid, a==b)

scheme = create.scheme("h1", spc, model,       parameters=       {"fem.solver.newton.linabstol": 1e-10,
        "fem.solver.newton.linreduction": 1e-10,
        "fem.solver.newton.verbose": 1,
        "fem.solver.newton.linear.verbose": 0})
# create a discrete solution over this space - will be initialized with zero by default

uh = create.function("discrete", spc, name="solution")


# In the following we implement a simple Newton solver: given an initial guess $u^0$ (here taken to be zero) solve for $n\geq 0$:
# \begin{align*}
#    u^{n+1} = u^n - DS(u^n)(S(u^n)-g)
# \end{align*}
# Where $g$ is a discrete function containing the boundary values in the Dirichlet nodes and zero otherwise.
#
# Let's first use the solve method on the scheme directly:

# In[11]:

_,info = scheme.solve(target = uh)
print("size:", grid.size(0), "newton iterations:", int(info['iterations']))
plot(grid, uh)


# Instead of `scheme.solve` we now use the call operator on the `scheme` (to compute $S(u^n$) as  well as `scheme.assemble` to get a copy of the system matrix in form of a scipy sparse row matrix. Note that this method is only available if the `storage` in the space is set `eigen`.

# In[12]:

# Let's first clear the solution again
uh.clear()
# Need to auxiliary function
res = uh.copy()

# Note: the following does not produce a copy of the dof
# vectors, but keep in mind that
# after grid adaptation the resulting numpy array
# will be invalid since the shared dof vector will have moved
# during its resizing - use copy=True to avoid this problem at
# the cost of a copy
sol_coeff = np.array( uh, copy=False )
res_coeff = np.array( res, copy=False )
n = 0

while True:
    scheme(uh, res)
    absF = sqrt( np.dot(res_coeff,res_coeff) )
    print("iterations ("+str(n)+")",absF)
    if absF < 1e-10:
        break
    matrix = scheme.assemble(uh)
    sol_coeff -= scipy.sparse.linalg.spsolve(matrix, res_coeff)
    n += 1

plot(grid, uh)


# We cam redo the above computation but now use the Newton solver available in sympy:

# In[13]:

# let's first set the solution back to zero - since it already contains the right values
uh.clear()
def f(x_coeff):
    x = spc.numpyfunction(x_coeff, "tmp")
    scheme(x,res)
    return res_coeff
# class for the derivative DS of S
class Df(scipy.sparse.linalg.LinearOperator):
    def __init__(self,x_coeff):
        self.shape = (sol_coeff.shape[0],sol_coeff.shape[0])
        self.dtype = sol_coeff.dtype
        # the following converts a given numpy array
        # into a discrete function over the given space
        x = spc.numpyfunction(x_coeff, "tmp")
        # store the assembled matrix
        self.jac = scheme.assemble(x)
    # reassemble the matrix DF(u) gmiven a dof vector for u
    def update(self,x_coeff,f):
        x = spc.numpyfunction(x_coeff, "tmp")
        # Note: the following does produce a copy of the matrix
        # and each call here will reproduce the full matrix
        # structure - no reuse possible in this version
        self.jac = scheme.assemble(x)
    # compute DS(u)^{-1}x for a given dof vector x
    def _matvec(self,x_coeff):
        return scipy.sparse.linalg.spsolve(self.jac, x_coeff)

# call the newton krylov solver from scipy
sol_coeff[:] = scipy.optimize.newton_krylov(f, sol_coeff,
            verbose=1, f_tol=1e-8,
            inner_M=Df(sol_coeff))

plot(grid, uh)


# ## Dirichlet boundary conditions
# demonstrate how to
# - set dirichlet={} in the model
# - set boundary ids in the grid construction

# ## DG Scheme
# show the
# - DG
