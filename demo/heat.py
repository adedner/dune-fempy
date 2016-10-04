from __future__ import print_function

import math
from ufl import *

from dune.grid import cartesianDomain
from dune.ufl import GridCoefficient, Space

import dune.create as create

# Crank Nicholson
theta = 0.5
deltaT = 0.01

def compute():
    # set up a 2d simplex grid over the interval [0,1]^2 with h = 1/16
    grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
    # set up a lagrange scalar space with polynomial order 2 over that grid
    spc = create.space("Lagrange", grid, dimrange=1, order=2)

    # set up initial conditions
    initial = lambda x: [ math.atan( (10.*x[0]*(1-x[0])*x[1]*(1-x[1]))**2 ) ]
    solution = spc.interpolate(initial, name="u")
    grid.writeVTK("heat", pointdata=[solution], number=0)

    # get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
    old_solution = spc.interpolate(initial, name="u_n")

    # now define the actual pde to solve:
    #            u - u_n deltaT laplace( theta u + (1-theta) u_n ) = 0
    uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    u_n = GridCoefficient(solution)
    tau = Constant(triangle)
    a = (inner(u - u_n, v) + tau * inner(grad(theta*u + (1-theta)*u_n), grad(v))) * dx

    # now generate the model code and compile
    model = create.model("elliptic", grid, a == 0, coefficients={u_n:old_solution})

    # setup structure for olver parameters
    solverParameter={"fem.solver.newton.linabstol": 1e-10,
                     "fem.solver.newton.linreduction": 1e-10,
                     "fem.solver.newton.verbose": 1,
                     "fem.solver.newton.linear.verbose": 1}
    # create the solver using a standard fem scheme
    scheme = create.scheme("h1", spc, model, "scheme", parameters=solverParameter)

    # now loop through time and output the solution after each time step
    steps = int(1 / deltaT)
    for n in range(1,steps+1):
        model.setConstant(tau,[deltaT])
        old_solution.assign(solution)
        scheme.solve(target=solution)
        grid.writeVTK("heat", pointdata=[solution], number=n)

compute()
