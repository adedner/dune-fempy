from __future__ import print_function

import math
from ufl import *

from dune.grid import cartesianDomain

import dune.create as create

from dune.ufl import DirichletBC, Space


# Crank Nicholson
# theta = 1.
deltaT = 0.001
viscosity = 1.0e0
import dune.fem
dune.fem.parameter.append({"fem.verboserank": 0,
                           "istl.preconditioning.method": "ilu",
                           "istl.preconditioning.iterations": 1,
                           "istl.preconditioning.relaxation": 1.2})

def compute():
    # set up a 2d simplex grid over the interval [0,1]^2 with h = 1/16
    # grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
    grid = create.grid("ALUConform", "../../data/channelcontraction.dgf", dimgrid=2)
    grid.hierarchicalGrid.globalRefine(6)
    # set up a lagrange scalar space with polynomial order 2 over that grid
    spc = create.space("lagrange", grid, dimrange=2, order=2, storage="istl")
    # spcP = create.space("lagrange", grid, dimrange=1, order=2)
    spcP = create.space("lagrange", grid, dimrange=1, order=1, storage="istl")

    # set up initial conditions

    velocity = spc.interpolate( [0,0], name="velocity")
    pressure = spcP.interpolate( lambda x: [0], name = "pressure", storage = "Istl" )

    solution = velocity, pressure

    vtk = grid.sequencedVTK("coupledNavstokes", pointvector={"velocity":velocity})
    vtk()

    # get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
    # old_solution = solution.copy()
    old_velo = velocity.copy()

    # old_solution.name = "uOld"
    old_velo.name = "uOld"

    # now define the actual pde to solve:
    #            u - u_n deltaT laplace( theta u + (1-theta) u_n ) = 0
    uflSpace = Space((grid.dimGrid, grid.dimWorld), 2)
    uflSpaceP = Space((grid.dimGrid, grid.dimWorld), 2)

    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    p = TrialFunction(uflSpaceP)
    q = TestFunction(uflSpaceP)
    x = SpatialCoordinate(uflSpace.cell())
    n = FacetNormal(uflSpace.cell())
    mu = 7.5 * 16

    dirichlet = conditional(abs(x[0])<1e-8,1.,0.)

    u_n = old_velo
    tau = Constant(triangle)
    a = (inner(u - u_n, v) + tau * viscosity * inner(grad(u), grad(v))) * dx
    a += -1*tau*inner(p[0],div(v)) * dx
    a += tau* inner(grad(u), outer(v, u)) * dx
    a += inner(div(u), q[0]) * dx


    # a -= (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
    # a += mu * inner(jump(u), jump(v)) * dS
    # a -= (inner(outer(u, n), grad(v)) + inner(grad(u), outer(v, n))) * ds
    # a += mu * inner(u, v) * ds
    # now generate the model code and compile
    #model = create.model("elliptic", grid, a == 0, coefficients={u_n:old_solution})
    #model.setConstant(tau,[deltaT])

    # model = create.model("integrands", grid, a == 0)
    model = create.model("elliptic", grid, a == 0,
            DirichletBC(uflSpace, [None,None], 2),       # bottom
            DirichletBC(uflSpace, [None,None], 3),      # top
            DirichletBC(uflSpace, [1,0], 4),          # left
            DirichletBC(uflSpace, [0,0], 1))                # right

    model.setConstant(tau, deltaT)

    # setup structure for olver parameters
    solverParameter={"fem.solver.newton.linabstol": 1e-13,
                     "fem.solver.newton.linreduction": 1e-13,
                     "fem.solver.newton.tolerance": 1e-12,
                     "fem.solver.newton.verbose": "true",
                     "fem.solver.newton.linear.verbose": "false"}
    # create the solver using a standard fem scheme
    # scheme = create.scheme("h1", spc, model, parameters=solverParameter)
    # scheme = create.scheme("h1galerkin", spc, model, parameters=solverParameter)
    # scheme = create.scheme("dggalerkin", spc, model, 15*theta*deltaT, parameters=solverParameter)

    scheme = create.scheme("h1", model, spc, parameters=solverParameter)

    # scheme = create.scheme("linearized", scheme, parameters=solverParameter)
    # scheme = create.scheme("linearized", scheme="h1", ubar=solution, space=spc, model=model, parameters=solverParameter)

    endTime = 0.05
    timeStep = deltaT
    time = timeStep
    while time < endTime:
        print( "Time is:", time )
        old_velo.assign(velocity)

        scheme.solve(target=velocity)
        vtk()
        time += timeStep

compute()
