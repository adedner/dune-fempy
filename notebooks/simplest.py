from __future__ import print_function, division

from ufl import as_vector, inner, grad, cos, pi, dx
from math import sqrt
from dune.ufl import DirichletBC

from dune.grid import structuredGrid
from dune.fem.space import lagrange
from dune.fem.plotting import plotPointData as plot
from dune.fem.function import integrate
from dune.fem.scheme import h1

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

# some postprocessing
plot(solution)
exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )
error = solution - exact
print("L^2 error:", sqrt( integrate(grid,error**2,order=5)[0] ) )
print("H^1 error:", sqrt( integrate(grid,inner(grad(error),grad(error)),order=5)[0] ) )

# heat equation
solution.clear()
un  = space.interpolate(-f,name="oldSolution")
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
