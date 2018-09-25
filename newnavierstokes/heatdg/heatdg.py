from __future__ import print_function

import math
from ufl import *

from dune.grid import cartesianDomain
from dune.ufl import Space

import dune.create as create

# Crank Nicholson
theta = 0.5
deltaT = 0.01

import dune.fem
dune.fem.parameter.append({"fem.verboserank": 0,
                           "istl.preconditioning.method": "ilu",
                           "istl.preconditioning.iterations": 1,
                           "istl.preconditioning.relaxation": 1.2})

def compute():
    # set up a 2d simplex grid over the interval [0,1]^2 with h = 1/16
    grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
    # set up a lagrange scalar space with polynomial order 2 over that grid
    spc = create.space("dgonb", grid, dimrange=1, order=2, storage="istl")
    # spc = create.space("dgonb", grid, dimrange=1, order=2)

    # set up initial conditions
    solution = spc.interpolate(lambda x: [math.atan((10.0 * x[0] * (1-x[0]) * x[1] * (1-x[1]))**2)], name="u")
    vtk = grid.sequencedVTK("heat-dg", pointdata=[solution])
    vtk()

    # get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
    old_solution = solution.copy();
    old_solution.name = "uOld"

    # now define the actual pde to solve:
    #            u - u_n deltaT laplace( theta u + (1-theta) u_n ) = 0
    uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    x = SpatialCoordinate(uflSpace.cell())
    n = FacetNormal(uflSpace.cell())
    mu = 7.5 * 16
    u_n = old_solution
    tau = Constant(triangle)
    a = (inner(u - u_n, v) + tau * inner(grad(theta*u + (1-theta)*u_n), grad(v))) * dx

    a -= (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
    a += mu * inner(jump(u), jump(v)) * dS
    a -= (inner(outer(u, n), grad(v)) + inner(grad(u), outer(v, n))) * ds
    a += mu * inner(u, v) * ds
    # now generate the model code and compile
    #model = create.model("elliptic", grid, a == 0, coefficients={u_n:old_solution})
    #model.setConstant(tau,[deltaT])

    model = create.model("integrands", grid, a == 0)
    model.setConstant(tau, deltaT)

    # setup structure for olver parameters
    solverParameter={"fem.solver.newton.linabstol": 1e-13,
                     "fem.solver.newton.linreduction": 1e-13,
                     "fem.solver.newton.tolerance": 1e-12,
                     "fem.solver.newton.verbose": "true",
                     "fem.solver.newton.linear.verbose": "false"}
    # create the solver using a standard fem scheme
    scheme = create.scheme("galerkin", spc, model, parameters=solverParameter)
    # scheme = create.scheme("h1galerkin", spc, model, parameters=solverParameter)
    # scheme = create.scheme("dggalerkin", spc, model, 15*theta*deltaT, parameters=solverParameter)

    # scheme = create.scheme("galerkin", model, spc, parameters=solverParameter)

    # scheme = create.scheme("linearized", scheme, parameters=solverParameter)
    # scheme = create.scheme("linearized", scheme="h1", ubar=solution, space=spc, model=model, parameters=solverParameter)

    # now loop through time and output the solution after each time step
    steps = int(0.1 / deltaT)
    for n in range(1,steps+1):
        old_solution.assign(solution)
        scheme.solve(target=solution)
        vtk()

compute()