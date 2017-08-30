# coding: utf-8

# # Different storage backends - using numpy/scipy [(Notebook)][1]
#
# [1]: _downloads/laplace-la.ipynb
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

# In[1]:

try:
    get_ipython().magic(u'matplotlib inline # can also use notebook or nbagg')
except:
    pass
from dune.generator import builder
import math
import numpy as np
import scipy.sparse.linalg
import scipy.optimize
import dune.grid
import dune.fem
from dune.fem.plotting import plotPointData as plot
import dune.create as create

from dune.ufl import Space
from ufl import TestFunction, TrialFunction, SpatialCoordinate, ds, dx, inner, grad

grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8]), dimgrid=2)

try:
    spc = create.space("Lagrange", grid, dimrange=1, order=1, storage='eigen')
except builder.ConfigurationError:
    exit(1)

d = 0.001
p = 1.7

uflSpace = Space(spc)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

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

# In[2]:

uh,info = scheme.solve(target = uh)
print("size:", grid.size(0), "newton iterations:", int(info['iterations']))
plot(uh)


# Instead of `scheme.solve` we now use the call operator on the `scheme` (to compute $S(u^n$) as  well as `scheme.assemble` to get a copy of the system matrix in form of a scipy sparse row matrix. Note that this method is only available if the `storage` in the space is set `eigen`.

# In[3]:

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
sol_coeff = uh.as_numpy
res_coeff = res.as_numpy
n = 0

while True:
    scheme(uh, res)
    absF = math.sqrt( np.dot(res_coeff,res_coeff) )
    print("iterations ("+str(n)+")",absF)
    if absF < 1e-10:
        break
    matrix = scheme.assemble(uh)
    sol_coeff -= scipy.sparse.linalg.spsolve(matrix, res_coeff)
    n += 1

plot(uh)


# We cam redo the above computation but now use the Newton solver available in sympy:

# In[7]:

# let's first set the solution back to zero - since it already contains the right values
uh.clear()
def f(x_coeff):
    x = spc.numpyFunction(x_coeff, "tmp")
    scheme(x,res)
    return res_coeff
# class for the derivative DS of S
class Df(scipy.sparse.linalg.LinearOperator):
    def __init__(self,x_coeff):
        self.shape = (sol_coeff.shape[0],sol_coeff.shape[0])
        self.dtype = sol_coeff.dtype
        # the following converts a given numpy array
        # into a discrete function over the given space
        x = spc.numpyFunction(x_coeff, "tmp")
        # store the assembled matrix
        self.jac = scheme.assemble(x)
    # reassemble the matrix DF(u) gmiven a dof vector for u
    def update(self,x_coeff,f):
        x = spc.numpyFunction(x_coeff, "tmp")
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

plot(uh)
