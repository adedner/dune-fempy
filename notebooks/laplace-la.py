
# coding: utf-8

# # Different storage backends - using numpy/scipy and petsc4py [(Notebook)][1]
#
# [1]: _downloads/laplace-la.ipynb
#
# For the first example we used solvers available in dune-fem - simple Krylov solvers with only diagonal preconditioning. Changing the `storage` argument in the construction of the space makes it possible to use more sophisticated solvers (either better preconditioners or direct solvers). For example
# ~~~
# spc = create.space("lagrange", grid, dimrange=1, order=1, storage="istl")
# ~~~
# in the above code will switch to the solvers from `dune-istl`, other options are for example `petsc` or `eigen`.
#
# Using the internal `fem` storage structure or the `eigen` matrix/vector strorage
# it is also possible to directly treate them as`numpy` vectors and an assembled system matrix can be stored in a `sympy` sparse matrix.
#
# __Note__: to use `eigen` matrices the `Eigen` package must be available and `dune-py` must have been configured with `Eigen`.
#
# Since we will be implementing a Newton solver first, let's study a truelly non linear problem - a version of the p-Laplace problem:
# \begin{gather}
#   - \frac{d}{2}\nabla\cdot |\nabla u|^{p-2}\nabla u + u = f
# \end{gather}

# In[ ]:


try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass
import sys
import petsc4py
from petsc4py import PETSc
petsc4py.init(sys.argv)

from dune.generator import builder
import math
import numpy as np
import scipy.sparse.linalg
import scipy.optimize
import dune.grid
from dune.fem.operator import linear
from dune.fem.plotting import plotPointData as plot
import dune.create as create

from dune.ufl import Space
from ufl import TestFunction, TrialFunction, SpatialCoordinate, ds, dx, inner, grad

grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8]), dimgrid=2)

spc = create.space("lagrange", grid, dimrange=1, order=1, storage='fem')

d = 0.001
p = 1.7

u = TrialFunction(spc)
v = TestFunction(spc)
x = SpatialCoordinate(spc.cell())

rhs = (x[0] + x[1]) * v[0]
a = (pow(d + inner(grad(u), grad(u)), (p-2)/2)*inner(grad(u), grad(v)) + inner(u, v)) * dx + 10*inner(u, v) * ds
b = rhs * dx + 10*rhs * ds
scheme = create.scheme("galerkin", a==b, spc, parameters=
       {"newton.tolerance": 1e-5, "newton.verbose": "false",
        "newton.linear.absolutetol": 1e-8, "newton.linear.reductiontol": 1e-8,
        "newton.linear.preconditioning.method": "ilu",
        "newton.linear.preconditioning.iterations": 1, "newton.linear.preconditioning.relaxation": 1.2,
        "newton.linear.verbose": "false"})
# create a discrete solution over this space - will be initialized with zero by default

uh = create.function("discrete", spc, name="solution")


# In the following we implement a simple Newton solver: given an initial guess $u^0$ (here taken to be zero) solve for $n\geq 0$:
# \begin{align*}
#    u^{n+1} = u^n - DS(u^n)(S(u^n)-g)
# \end{align*}
# Where $g$ is a discrete function containing the boundary values in the Dirichlet nodes and zero otherwise.
#
# Let's first use the solve method on the scheme directly:

# In[ ]:


info = scheme.solve(target = uh)
print(info)
plot(uh)


# Instead of `scheme.solve` we now use the call operator on the `scheme` (to compute $S(u^n$) as  well as `scheme.assemble` to get a copy of the system matrix in form of a scipy sparse row matrix. Note that this method is only available if the `storage` in the space is set `eigen`.

# In[ ]:


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

matrix = linear(scheme)
matrix_coeff = matrix.as_numpy

import pyamg
while True:
    scheme(uh, res)
    absF = math.sqrt( np.dot(res_coeff,res_coeff) )
    print("iterations ("+str(n)+")",absF)
    if absF < 1e-10:
        break
    scheme.jacobian(uh,matrix)
    ### spsolve fails because it changes the numpy matrix (permutes the col)
    ### and that fails for some reason (not yet clear why)
    # matrix_coeff = matrix.as_numpy
    # sol_coeff -= scipy.sparse.linalg.spsolve(matrix_coeff_tmp, res_coeff)
    sol_coeff -= scipy.sparse.linalg.cg(matrix_coeff, res_coeff, tol=1e-10)[0]

    # print("setting up amg")
    # construct the multigrid hierarchy
    # ml = pyamg.ruge_stuben_solver(matrix_coeff)
    # ml = pyamg.smoothed_aggregation_solver(matrix_coeff)
    # ml = pyamg.rootnode_solver(matrix_coeff, smooth=('energy', {'degree':2}), strength='evolution' )
    # print(ml)
    # M = ml.aspreconditioner(cycle='V')
    # print("solving...")
    # sol_coeff -= scipy.sparse.linalg.cg(matrix_coeff, res_coeff, tol=1e-10, M=M)[0]
    # print("...done")

    # print("setting up amg")
    # A = scipy.sparse.csr.csr_matrix(matrix_coeff.todense())
    # ml = pyamg.smoothed_aggregation_solver(A, max_coarse=10)
    # print(ml)
    # residuals = []
    # print("solving...")
    # sol_coeff -= ml.solve(res_coeff, tol=1e-10, accel='cg', residuals=residuals)
    # residuals = residuals / residuals[0]
    # print("...done")

    n += 1

plot(uh)


# We cam redo the above computation but now use the Newton solver available in sympy:

# In[ ]:


# let's first set the solution back to zero - since it already contains the right values
uh.clear()
def f(x_coeff):
    x = spc.function("tmp", dofVector=x_coeff)
    scheme(x,res)
    return res_coeff
# class for the derivative DS of S
class Df(scipy.sparse.linalg.LinearOperator):
    def __init__(self,x_coeff):
        self.shape = (sol_coeff.shape[0],sol_coeff.shape[0])
        self.dtype = sol_coeff.dtype
        # setup the assembled matrix
        self.linOp = linear(scheme)
    # reassemble the matrix DF(u) gmiven a dof vector for u
    def update(self,x_coeff,f):
        # the following converts a given numpy array
        # into a discrete function over the given space
        x = spc.function("tmp", dofVector=x_coeff)
        scheme.jacobian(x,self.linOp)
    # compute DS(u)^{-1}x for a given dof vector x
    def _matvec(self,x_coeff):
        return scipy.sparse.linalg.cg(self.linOp.as_numpy, x_coeff, tol=1e-10)[0]

# call the newton krylov solver from scipy
sol_coeff[:] = scipy.optimize.newton_krylov(f, sol_coeff,
            verbose=1, f_tol=1e-8,
            inner_M=Df(sol_coeff))

plot(uh)

# import sys
# import petsc4py
# from petsc4py import PETSc
# petsc4py.init(sys.argv)
# sys.exit(0)


# We can also use the package `petsc4py` to solve the problem.
#
# __Note__: make sure that `dune` has been configured using the same version of `petsc` used for `petsc4py`
#
# The first step is to change the storage in the space. Since also requires setting up the scheme and siscrete functions again to use the new storage structure.
#
# We can directly use the `petsc` solvers by invoking `solve` on the scheme as before.

# In[ ]:


try:
    # import petsc4py, sys
    # from petsc4py import PETSc
    # petsc4py.init(sys.argv)
    spc = create.space("lagrange", grid, dimrange=1, order=1, storage='petsc')
    scheme = create.scheme("galerkin", a==b, spc,
                            parameters={"petsc.preconditioning.method":"sor"})
    # first we will use the petsc solver available in the `dune-fem` package (using the sor preconditioner)
    uh   = spc.interpolate([0],name="petsc")
    info = scheme.solve(target=uh)
    print(info)
    plot(uh)
except ImportError:
    print("petsc4py could not be imported")
    petsc4py = False

# Next we will implement the Newton loop in Python using `petsc4py` to solve the linear systems
# Need to auxiliary function and set `uh` back to zero.
# We can access the `petsc` vectors by calling `as_petsc` on the discrete function. Note that this property will only be available if the discrete function is an element of a space with storage `petsc`.
# The method `assemble` on the scheme now returns the sparse `petsc` matrix and so we can direclty use the `ksp` class from `petsc4py`:

# In[ ]:


if petsc4py:
    uh.clear()
    res = uh.copy()

    sol_coeff = uh.as_petsc
    res_coeff = res.as_petsc

    ksp = PETSc.KSP()
    ksp.create(PETSc.COMM_WORLD)
    # use conjugate gradients method
    ksp.setType("cg")
    # and incomplete Cholesky
    ksp.getPC().setType("icc")

    n = 0
    linOp  = linear(scheme)
    ksp.setOperators(linOp.as_petsc)
    ksp.setFromOptions()
    while True:
        scheme(uh, res)
        absF = math.sqrt( res_coeff.dot(res_coeff) )
        print("iterations ("+str(n)+")",absF)
        if absF < 1e-10:
            break
        scheme.jacobian(uh,linOp)
        ksp.solve(res_coeff, res_coeff)
        sol_coeff -= res_coeff
        n += 1
    plot(uh)


# Finally we weill use `petsc`'s non-linear solvers (the `snes` classes) directly:

# In[ ]:


if petsc4py:
    uh.clear()
    res.clear()
    def f(snes, X, F):
        inDF  = spc.function("tmp", dofVector=X)
        outDF = spc.function("tmp", dofVector=F)
        scheme(inDF,outDF)
    def Df(snes, x, m, b):
        inDF = spc.function("tmp", dofVector=x)
        scheme.jacobian(inDF, linOp)
        return PETSc.Mat.Structure.SAME_NONZERO_PATTERN

    matrix = linOp.as_petsc
    snes = PETSc.SNES().create()
    snes.setMonitor(lambda snes,i,r: print(i,r,flush=True))
    b = res_coeff.duplicate()
    snes.setFunction(f, b)
    snes.setUseMF(False)
    snes.setJacobian(Df,matrix,matrix)
    snes.getKSP().setType("cg")
    snes.setFromOptions()
    snes.solve(res_coeff, sol_coeff)
    plot(uh)


# __Note__:
# The method `as_numpy, as_petsc` returning the `dof` vector either as a `numpy` or a `petsc` do not lead to a copy of the data and the same is true for the `function` method on the space. In the `numpy` case we can use `Python`'s buffer protocol to use the same underlying storage. In the case of `petsc` the underlying `Vec` can be shared. In the case of matrices the situation is not yet as clear: `scheme.assemble` returns a copy of the data in the `scipy` case while the `Mat` structure is shared between `c++` and  `Python` in the `petsc` case. But at the time of writting it is not possible to pass in the `Mat` structure to the `scheme.assemble` method from the outside. That is why it is necessary to copy the data when using the `snes` non linear solver as seen above.
