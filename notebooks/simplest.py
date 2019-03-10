from __future__ import print_function, division

import numpy as np
from ufl import as_vector, inner, grad, sin, cos, pi, dx, atan_2, conditional, dot
from math import sqrt, pi
from dune.ufl import DirichletBC, Constant

import dune.create as create
from dune.grid import structuredGrid
from dune.fem.space import lagrange
from dune.fem.plotting import plotPointData as plot
from dune.fem.function import integrate
from dune.fem.scheme import galerkin

import ufl

def plot(*args,**kwargs):
    pass

###################################################################

# grid and space
grid  = structuredGrid([0, 0], [1, 1], [16, 16])
space = lagrange(grid, order=1)

# ufl
u = ufl.TrialFunction(space)
v = ufl.TestFunction(space)
x = ufl.SpatialCoordinate(space.cell())

f = (8*pi*pi+1)*cos(2*pi*x[0])*cos(2*pi*x[1])

# elliptic equation
scheme = galerkin( (
         u*v  + dot(grad(u),grad(v)) )*dx == f*v*dx )

solution = space.interpolate([0],name="solution")
info = scheme.solve(target=solution)

# some postprocessing
plot(solution)
exact = cos(2.*pi*x[0])*cos(2.*pi*x[1])
error = solution - exact
print("L^2 and H^1 error:",
  [ sqrt(e) for e in
    integrate(grid,[error**2,inner(grad(error),grad(error))], order=5) ] )
plot(error,grid=grid)

# heat equation
dt = 0.01
solution.clear()
un  = space.interpolate(-exact,name="oldSolution")
tau = Constant(dt,name="tau")
scheme = galerkin( ( u*v  + tau*dot(grad(u),grad(v)) )*dx == (un+tau*f)*v*dx )

# compute until stationary solution reached (same final solution as above)
t  = 0
dist = integrate(grid,error**2,order=5)[0]
while dist > 0.01:
    info = scheme.solve(target=solution)
    t += dt
    un.assign(solution)
    dist = integrate(grid,error**2,order=5)[0]
    print(t,dist,info)
plot(solution)

# complicated domain with boundary conditions
vertices = np.array([(0,0), (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1)])
triangles = np.array([(0,1,2), (0,2,3), (0,3,4), (0,4,5), (0,5,6), (0,6,7)])
grid = create.grid("ALUConform", {"vertices": vertices, "simplices": triangles}, dimgrid=2)
grid.hierarchicalGrid.globalRefine(4)
space = lagrange(grid, dimRange=1, order=1)
u = ufl.TrialFunction(space)
v = ufl.TestFunction(space)
x = ufl.SpatialCoordinate(space.cell())
phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*pi, 0)
exact = as_vector([inner(x,x)**(0.5*180./270.) * sin((180./270.) * phi)])
a = inner(grad(u), grad(v))*dx
scheme = galerkin([a==0,DirichletBC(space,exact,1)])
solution = space.interpolate([0],name="solution")
scheme.solve(target=solution)
plot(solution)
error = solution - exact
print("L^2 error:", sqrt( integrate(grid,error**2,order=5)[0] ) )
# problem with division by zero in the following
print("H^1 error:", sqrt( integrate(grid,inner(grad(error),grad(error)),order=5)[0] ) )
