"""Solve the Laplace equation
"""
from __future__ import print_function

from mpi4py import MPI

import math,sys
from ufl import *
import scipy
import scipy.sparse.linalg
import scipy.optimize
import numpy

import dune.ufl
import dune.create as create
import dune.fem

# bug fix:
# sol.assign(bnd)  why does this not work - the below does
# sol_coeff[:] = bnd_coeff[:]
# not possible? sol = spc.interpolate("solution", [0])
# not possible? rhs = spc.interpolate("solution", [0])

def compute():
    dune.fem.parameter.append("../data/parameter")

    dgf = """
            Interval
            0   0                                   % first corner
            1   1                                   % second corner
            8   8                                   % 8 cells in each direction
            #
            BoundaryDomain
            default 1                               % all boundaries have id 1
            #
         """

    grid = create.grid("ALUConform", dune.grid.string2dgf(dgf), dimgrid=2)
    spc  = create.space("Lagrange", grid, dimrange=1, order=2)

    uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1, field="double")
    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    x = SpatialCoordinate(uflSpace.cell())

    exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )

    a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
    a = a + 20./(u[0]*u[0]+1.) * v[0] * dx

    model = create.model("elliptic", grid, a==0, exact=exact, dirichlet={ 1:exact } )

    # scheme = dune.fem.create.scheme("DGFemScheme", spc, model,\
    scheme = create.scheme("h1", spc, model,\
           parameters=\
           {"fem.solver.newton.tolerance": 1e-10,
            "fem.solver.newton.linabstol": 1e-12,
            "fem.solver.newton.linreduction": 1e-12,
            "fem.solver.newton.verbose": 1,
            "fem.solver.newton.linear.verbose": 1},\
            storage="eigen")

    exact_gf = create.function("ufl", grid, "exact", 5, exact)

    sol = create.discretefunction("eigen", spc, name="solution")
    sol.interpolate( [0] )
    for i in range(2):
        print("solve on level",i)
        uh = scheme.solve()

        def ver1():
            # use the apply/assemble methods directly in a simple Newton loop
            # - exact copy of simple version in dune-fem
            sol.clear()
            # Note: the following does not produce a copy of the dof
            # vectors, but keep in mind that
            # after grid adaptation the resulting numpy array
            # will be invalid since the shared dof vector will have moved
            # during its resizing - use copy=True to avoid this problem at
            # the cose of a copy
            sol_coeff = numpy.array( sol, copy=False )
            res = spc.interpolate(lambda _,x:[0], name="res", storage="eigen")
            res_coeff = numpy.array( res, copy=False )
            bnd = spc.interpolate(lambda _,x:[0], name="bnd", storage="eigen")
            bnd_coeff = numpy.array( bnd, copy=False )
            scheme.constraint(bnd)   # note: no rhs functional allowed
            n = 0
            while True:
                scheme(sol, res)
                res_coeff -= bnd_coeff
                absF = math.sqrt( numpy.dot(res_coeff,res_coeff) )
                print("Numpy iterations ("+str(n)+")",absF)
                if absF < 1e-10:
                    break

                matrix = scheme.assemble(sol)
                sol_coeff -= scipy.sparse.linalg.spsolve(matrix, res_coeff)
                n += 1

        def ver2():
            # use the numpy's newton_krylov method
            sol_coeff = numpy.array( sol, copy=False )
            sol.clear()
            rhs = spc.interpolate(lambda x:[0], name="rhs", storage="eigen")
            rhs_coeff = numpy.array( rhs, copy=False )
            bnd = spc.interpolate(lambda x:[0], name="bnd", storage="eigen")
            bnd_coeff = numpy.array( bnd, copy=False )
            scheme.constraint(bnd)
            def f(x_coeff):
                x = spc.numpyfunction(x_coeff, "tmp")
                scheme(x,rhs)
                return rhs_coeff - bnd_coeff
            class Df(scipy.sparse.linalg.LinearOperator):
                def __init__(self,x_coeff):
                    self.shape = (sol_coeff.shape[0],sol_coeff.shape[0])
                    self.dtype = sol_coeff.dtype
                    x = spc.numpyfunction(x_coeff, "tmp")
                    # x = create.function("numpy", spc, x_coeff, "tmp")
                    self.jac = scheme.assemble(x)
                def update(self,x_coeff,f):
                    x = spc.numpyfunction(x_coeff, "tmp")
                    # Note: the following does produce a copy of the matrix
                    # and each call here will reproduce the full matrix
                    # structure - no reuse possible in this version.
                    # Also: the assemble method is only available on
                    # schemes with storage="eigen"
                    self.jac = scheme.assemble(x)
                def _matvec(self,x_coeff):
                    return scipy.sparse.linalg.spsolve(self.jac, x_coeff)
            sol_coeff[:] = scipy.optimize.newton_krylov(f, sol_coeff,
                        verbose=1, f_tol=1e-10,
                        inner_M=Df(sol_coeff))

        ver1()
        ver2()

        # error between fem and exact solution
        def l2error(en,x):
            val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
            return [ val[0]*val[0] ];
        l2error_gf = create.function("local", grid, "error", 5, l2error )
        error = math.sqrt( l2error_gf.integrate()[0] )
        print("size:",grid.size(0),"L2-error:",error)
        grid.writeVTK("laplace", pointdata=[ uh,l2error_gf,sol ], number=i)

        # error between numpy.newton_krylov and exact solution
        def l2error(en,x):
            val = sol.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
            return [ val[0]*val[0] ];
        l2error_gf = create.function("local", grid, "error", 5, l2error )
        error = math.sqrt( l2error_gf.integrate()[0] )
        print("size:",grid.size(0),"L2-error:",error)

        grid.hierarchicalGrid.globalRefine(2)

compute()
