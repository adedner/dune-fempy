"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.

  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary

  u_D = 1 + x^2 + 2y^2
    f = -6
"""

from __future__ import print_function
from dune.fenics import *
import dune.fem
dune.fem.parameter.append({"fem.verboserank": 0})

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)

x = SpatialCoordinate( V )
# Define boundary condition
# original
# --------
# u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
# def boundary(x, on_boundary):
#     return on_boundary
# u_D = Expression(’1 + x[0]*x[0] + 2*x[1]*x[1]’, degree=2)
# bc = DirichletBC(V, u_D, boundary)
# dune-fempy:
# We don't have `Expression` i.e. generating code directly from a string
# Was implemented some time ago and could be revived (probably)
# What would work now is - which I think is actually better :-)
# What the Expression class allows you to do is
# u_D.t = t (see heat example page 41/42
# We could define an `Expression` class which looks similar, definig
# `NamedConstants` for the argument like `t` and return a uflFunction - it
# would still take a UFL expression instead of a string
x = SpatialCoordinate(mesh)
u_D = 1 + x[0]*x[0] + 2*x[1]*x[1]
bc = [DirichletBC(V, u_D, i) for i in range(1,5)] # boundary ids are 1,..,4 or yasp

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

f = Constant(-6.0)
a = dot(grad(u), grad(v)) * dx
L = f*v*dx

# Compute solution
u = Function(V)
# default solver in fenics is LU decompose (umfpack/petsc.lu or something like that?)
# we could try testing for existence of umpack/petsc and use that but it should
# not be 'cg' by default if everything else is not available
solve(a == L, u, bc)

# Plot solution and mesh
plot(u)
plot(mesh)

# Save solution to file in VTK format
# NOTE: might be possible - need to define a `File` class and add an
# method __lshift__(self,other) - not quite sure about the semantics, is
# the file written when executing `vtkfile << u`? Or is the pvd (an xml
# file format) simply appended to? Could be quite nice to have but that
# `File` class is perhaps best written in python...
# vtkfile = File('poisson/solution.pvd')
# vtkfile << u

# Compute error in L2 norm
error_L2 = errornorm(u_D, u, 'L2')

# Compute maximum error at vertices
# original
# --------
# vertex_values_u_D = u_D.compute_vertex_values(mesh)
# vertex_values_u = u.compute_vertex_values(mesh)
# import numpy as np
# error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
# fempy
# -----
# we can't do the first one since u_D in our case is simply a UFL
# expression without any change - so need something like `Expression` but
# without the string
# vertex_values_u_D = u_D.compute_vertex_values.pointData()
vertex_values_u = u.pointData()
import numpy as np
error_max = np.max(np.abs(vertex_values_u))

# Print errors
print('error_L2  =', error_L2)
print('error_max =', error_max)

# Hold plot
interactive()
