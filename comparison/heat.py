from __future__ import print_function

import imp
import math

from dune.grid import cartesianDomain

import dune.create as create

import dune.fem
from dune.fem import parameter
parameter.append( "parameter" )
dt = 0.02

def compute():
    # set up a 2d simplex grid over the interval [0,1]^2 with h = 1/16
    grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[120,120]), dimgrid=2)
    # set up a lagrange scalar space with polynomial order 2 over that grid
    spc = create.space("lagrange", grid, dimrange=1, order=2, storage="fem")

    # set up initial conditions
    solution = spc.interpolate(lambda x: [math.atan((10.0 * x[0] * (1-x[0]) * x[1] * (1-x[1]))**2)], name="u")
    grid.writeVTK("heat", pointdata=[solution], number=0)

    # get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
    u_n = solution.copy();

    heatUFL = imp.load_source("heatUFL", "./heat.ufl")

    model = create.model("elliptic", grid, heatUFL.F==0,
                         coefficients={heatUFL.u_n:u_n})
    model.tau = dt
    model.theta = 0.5

    # create the solver using a standard fem scheme
    scheme = create.scheme("h1", spc, model, solver="cg")

    # now loop through time and output the solution after each time step
    steps = int(1. / dt)
    for n in range(0,steps):
        model.time = n*dt
        u_n.assign(solution)
        scheme.solve(target=solution)

    grid.writeVTK("heat", pointdata=[solution], number=n)

compute()
