from __future__ import print_function

from mpi4py import MPI

import math
from ufl import *

import dune.models.elliptic
import dune.ufl
import dune.fem

# Crank Nicholson
theta = 0.5
deltaT = 0.01

# set up a 2d simplex grid over the interval [0,1]^2 with h = 1/16
grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
# set up a lagrange scalar space with polynomial order 2 over that grid
spc = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

# now define the actual pde to solve:
#            u - u_n deltaT laplace( theta u + (1-theta) u_n ) = 0
uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
u_n = Coefficient(uflSpace)
a = (inner(u - u_n, v) + deltaT * inner(grad(theta*u + (1-theta)*u_n), grad(v))) * dx

# now generate the model code and compile
model = dune.models.elliptic.importModel(grid, dune.models.elliptic.compileUFL(a == 0)).get()

# create the solver using a standard fem scheme
scheme = dune.fem.create.scheme("FemScheme", spc, model, "scheme")

# set up initial conditions
solution = spc.interpolate(lambda x: [ math.atan( (10.*x[0]*(1-x[0])*x[1]*(1-x[1]))**2 ) ], name="u")
grid.writeVTK("heat", pointdata=[solution], number=0)

# get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
old_solution = spc.interpolate(lambda x: [ math.atan( (10.*x[0]*(1-x[0])*x[1]*(1-x[1]))**2 ) ], name="u_n")
model.setCoefficient(u_n, old_solution)

# now loop through time and output the solution after each time step
steps = int(1 / deltaT)
for n in range(1,steps+1):
    old_solution.assign(solution)
    scheme.solve(target=solution)
    grid.writeVTK("heat", pointdata=[solution], number=n)
