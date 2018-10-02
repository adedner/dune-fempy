from __future__ import print_function

import math
from ufl import *

from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, dot, div, grad, dx, as_vector, transpose, Identity

from dune.grid import cartesianDomain

from dune.ufl import NamedConstant

import dune.create as create

from dune.ufl import DirichletBC, Space


# Crank Nicholson
theta = 1
deltaT = 0.03
viscosity = 1.0e-1
import dune.fem
dune.fem.parameter.append({"fem.verboserank": 0,
                           "istl.preconditioning.method": "ilu",
                           "istl.preconditioning.iterations": 1,
                           "istl.preconditioning.relaxation": 1.2})

def compute():
    # set up a 2d simplex grid over the interval [0,1]^2 with h = 1/16
    # grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
    grid = create.grid("ALUConform", "../../data/channelcontraction.dgf", dimgrid=2)
    grid.hierarchicalGrid.globalRefine(9)
    # set up a lagrange scalar space with polynomial order 2 over that grid
    spc = create.space("dgonb", grid, dimrange=2, order=2, storage="istl")
    # spcP = create.space("lagrange", grid, dimrange=1, order=2)
    spcP = create.space("dgonb", grid, dimrange=1, order=1, storage="istl")

    # set up initial conditions

    velocity = spc.interpolate( [0,0], name="velocity")
    pressure = spcP.interpolate( lambda x: [0], name = "pressure", storage = "Istl" )

    solution = velocity, pressure

    vtk = grid.sequencedVTK("coupledNavstokes", pointdata={"velocity":velocity})
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
    # mu = 7.5 * 16
    mu = 10e4


    phi = conditional(x[0] < 1e-8, 1, 0)

    u_d = as_vector( [phi, 0] )


    dirichlet = conditional(x[0]<1e-8,1.,0. ) + conditional(x[0]>29.99 ,1.,0. )

    # u_n = old_velo
    # tau = Constant(triangle)
    # a = (inner(u - u_n, v) + tau * viscosity * inner(grad(theta*u + (1-theta)*u_n), grad(v))) * dx
    # a += -1*tau*inner(p[0],div(v)) * dx
    # a += tau* inner(grad(u), outer(v, u)) * dx
    # a += div(u)*q[0] * dx


    tau = NamedConstant(spc,"tau")

    # diffusion
    a  = viscosity * inner(grad(u), grad(v)) * dx
    a -= tau * viscosity * (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
    a += tau * mu *  viscosity * inner(jump(u), jump(v)) * dS
    # a -= tau * (inner(outer(u_d, n), grad(v)) + inner(grad(u_d), outer(v, n))) * ds
    a -= tau *  viscosity * (inner(outer(u_d, n), grad(v)) ) * dirichlet *  ds
    a += tau * viscosity * mu * inner(u_d, v) * dirichlet * ds

    # advection
    a += tau* inner(grad(u), outer(v, u)) * dx
    a += 0.5 * tau* inner(v, outer(div(u), u)) * dx
    a -= tau * (inner(outer(jump(u), n('+')), outer(avg(u),avg(v)) )) * dS
    # a -= 0.5 * tau * (inner(outer(jump(u), n('+')), avg(u*v) )) * dS

    # pressure part
    a += -1*tau*inner(p[0],div(v)) * dx
    a += tau * (inner(outer(avg(p[0]), n('+')), jump(v))) * dS
    # a += tau * inner(p_d, v) * dirichlet * ds

    # incompressibilty (pressure part)
    a -= div(u)*q[0] * dx
    a += (inner(outer(avg(q[0]), n('+')), jump(u))) * dS
    a += (inner(outer(q[0], n), u_d)) * dirichlet * ds
    # a += (1/viscosity) * mu * inner(jump(u), outer(jump(q[0]), n('+'))) * dS
    # a += (1/viscosity) * mu * inner(u_d, outer(q[0], n)) *dirichlet *  ds

    # complete coupled scheme
    a = inner(u - old_velo, v) * dx + tau*theta*a + tau*(1-theta)*replace(a,{u:old_velo})


    # a -= (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
    # a += mu * inner(jump(u), jump(v)) * dS
    # a -= (inner(outer(u, n), grad(v)) + inner(grad(u), outer(v, n))) * ds
    # a += mu * inner(u, v) * ds
    # now generate the model code and compile
    #model = create.model("elliptic", grid, a == 0, coefficients={u_n:old_solution})
    #model.setConstant(tau,[deltaT])

    model = create.model("integrands", grid, a == 0)
    # model = create.model("elliptic", grid, a == 0,
    #         DirichletBC(uflSpace, [None,None], 2),       # bottom
    #         DirichletBC(uflSpace, [None,None], 3),      # top
    #         DirichletBC(uflSpace, [1,0], 4),          # left
    #         DirichletBC(uflSpace, [0,0], 1))                # right

    model.setConstant(tau, deltaT)

    # setup structure for olver parameters
    solverParameter={"fem.solver.newton.linabstol": 1e-12,
                     "fem.solver.newton.linreduction": 1e-12,
                     "fem.solver.newton.tolerance": 1e-11,
                     "fem.solver.newton.verbose": "true",
                     "fem.solver.newton.linear.verbose": "false"}
    # create the solver using a standard fem scheme
    # scheme = create.scheme("h1", spc, model, parameters=solverParameter)
    # scheme = create.scheme("h1galerkin", spc, model, parameters=solverParameter)
    # scheme = create.scheme("dggalerkin", spc, model, 15*theta*deltaT, parameters=solverParameter)

    scheme = create.scheme("galerkin", model, spc, parameters=solverParameter)

    # scheme = create.scheme("linearized", scheme, parameters=solverParameter)
    # scheme = create.scheme("linearized", scheme="h1", ubar=solution, space=spc, model=model, parameters=solverParameter)

    # now loop through time and output the solution after each time step
    steps = 5000
    for n in range(1,steps+1):
        # old_solution.assign(solution)
        old_velo.assign(velocity)

        scheme.solve(target=velocity)
        vtk()

compute()
