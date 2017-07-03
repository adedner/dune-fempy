from __future__ import print_function, division

import numpy as np
from ufl import as_vector, inner, grad, sin, cos, pi, dx, atan_2, conditional
from math import sqrt, pi
from dune.ufl import DirichletBC

import dune.create as create
from dune.grid import structuredGrid
from dune.fem.space import lagrange
from dune.fem.plotting import plotPointData as plot
from dune.fem.function import integrate
from dune.fem.scheme import h1

import ufl

def plot(*args,**kwargs):
    pass

###################################################################

# grid and space
grid  = structuredGrid([0, 0], [1, 1], [16, 16])
space = lagrange(grid, dimrange=1, order=1)

# ufl
u = space.uflTrialFunction()
v = space.uflTestFunction()
x = space.uflSpatialCoordinate()
f = as_vector( [(8*pi*pi+1)*cos(2*pi*x[0])*cos(2*pi*x[1])] )

# elliptic equation
scheme = h1(space,
          ( inner(u,v)  + inner(grad(u),grad(v)) )*dx == inner(f,v)*dx )

solution, info = scheme.solve()

print(space.size)
uh_dofsA = solution.as_numpy
uh_dofsB = solution.dofVector
uh_dofsC = np.array( solution.dofVector, copy=False )
print(solution.array.size,solution.array.strides)
print(uh_dofsC.size,uh_dofsC.strides)
print(uh_dofsB.size)
for a in range(space.size):
    print(a,uh_dofsB[a])
    print(a,uh_dofsC[a])
    print(a,uh_dofsA[a])
    print(a,uh_dofsA[a]-uh_dofsB[a],uh_dofsA[a]-uh_dofsC[a])
for a in range(space.size):
    uh_dofsB[a] = 0
for a in range(space.size):
    print(a,uh_dofsA[a],uh_dofsB[a],uh_dofsC[a])
for a in range(space.size):
    uh_dofsA[a] = 1
for a in range(space.size):
    print(a,uh_dofsA[a],uh_dofsB[a],uh_dofsC[a])
print(uh_dofsA.size, uh_dofsA.strides)
print(len(uh_dofsB),uh_dofsB.size)
print(uh_dofsC.size, uh_dofsC.strides)
exit(1)

print(isinstance(solution,ufl.Coefficient))
print(dir(solution))
print(type(solution))
print(solution.space)
help(solution)
# help(solution.GridFunction)
# help(solution.Coefficient)
# help(solution.grid)
# help(solution.interpolate)

# some postprocessing
plot(solution)
exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )
error = solution - exact
print("L^2 error:", sqrt( integrate(grid,error**2,order=5)[0] ) )
print("H^1 error:", sqrt( integrate(grid,inner(grad(error),grad(error)),order=5)[0] ) )


# heat equation
solution.clear()
un  = space.interpolate(-exact,name="oldSolution")
tau = space.uflNamedConstant(name="tau")
scheme = h1(space,
        ( inner(u,v)  + tau*inner(grad(u),grad(v)) )*dx == inner(un+tau*f,v)*dx )

# compute until stationary solution reached (same final solution as above)
t  = 0
dt = 0.01
scheme.model.tau = dt
dist = integrate(grid,error**2,order=5)[0]
while dist > 0.01:
    info = scheme.solve(target=solution)
    t += dt
    un.assign(solution)
    dist = integrate(grid,error**2,order=5)[0]
    print(t,dist)
plot(solution)

# complicated domain with boundary conditions
vertices = np.array([(0,0), (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1)])
triangles = np.array([(0,1,2), (0,2,3), (0,3,4), (0,4,5), (0,5,6), (0,6,7)])
grid = create.grid("ALUConform", {"vertices": vertices, "simplices": triangles}, dimgrid=2)
grid.hierarchicalGrid.globalRefine(4)
space = lagrange(grid, dimrange=1, order=1)
u = space.uflTrialFunction()
v = space.uflTestFunction()
x = space.uflSpatialCoordinate()
phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*pi, 0)
exact = as_vector([inner(x,x)**(0.5*180./270.) * sin((180./270.) * phi)])
a = inner(grad(u), grad(v))*dx
scheme = h1(space, [a==0,DirichletBC(space.uflSpace(),exact,1)])
solution, _ = scheme.solve()
plot(solution)
error = solution - exact
print("L^2 error:", sqrt( integrate(grid,error**2,order=5)[0] ) )
# problem with division by zero in the following
# print("H^1 error:", sqrt( integrate(grid,inner(grad(error),grad(error)),order=5)[0] ) )
