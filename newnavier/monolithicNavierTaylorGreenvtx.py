from __future__ import print_function

import math
from ufl import *


from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, dot, div, grad, dx, as_vector, transpose, Identity

from dune.grid import cartesianDomain, Marker

import dune.create as create

from dune.ufl import DirichletBC, NamedConstant

from ufl import as_vector, dx, grad, inner, replace, exp, dot
from dune.fem.function import levelFunction

from dune.fem import parameter, adapt

import dune.fem as fem
from math import pi,log10
from dune.fem.function import integrate

# Crank Nicholson
# theta = 1.
deltaT = 1e-3 # 0.001
viscosity = 1.0e-2
d = 0.00
pnb = 2.0
maxLevel = 8

import dune.fem
dune.fem.parameter.append({"fem.verboserank": 0,
                           "istl.preconditioning.method": "ilu",
                           "istl.preconditioning.iterations": 1,
                           "istl.preconditioning.relaxation": 1.2})

def compute():
    # set up a 2d simplex grid over the interval [0,1]^2 with h = 1/16
    # grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
    # grid = create.grid("ALUCube", "../../data/channelexpansion.dgf", dimgrid=2)
    grid = create.grid("ALUCube",constructor=cartesianDomain([-1,-1],[1,1],[50,50]))
    # grid.hierarchicalGrid.globalRefine(6)
    # set up a lagrange scalar space with polynomial order 2 over that grid
    spc = create.space("lagrange", grid, dimrange=3, order=2, storage="istl")

    # set up initial conditions


    # now define the actual pde to solve:
    #            u - u_n deltaT laplace( theta u + (1-theta) u_n ) = 0
    cell  = spc.cell()

    trial = TrialFunction(spc)
    test  = TestFunction(spc)
    u = as_vector([trial[0],trial[1]])
    v = as_vector([test[0],test[1]])
    p = as_vector([trial[2]])
    q = as_vector([test[2]])
    x = SpatialCoordinate(spc)

    # mu    = NamedConstant(cell, "mu")
    # nu    = NamedConstant(cell, "nu")
    time     = NamedConstant(triangle, "t")     # current time

    mu = 0.01
    nu = 0.0

    gamma_t = exp(-2.*pi*pi*mu*time)
    exact_u     = as_vector( [-1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t, sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*mu*time),  -0.25*(cos(2.*pi*x[0])+cos(2.*pi*x[1]))*exp(-4.*pi*pi*mu*time)])
    exact_p     = as_vector( [ -0.25*(cos(2.*pi*x[0])+cos(2.*pi*x[1]))*exp(-4.*pi*pi*mu*time) ] )

    exact_uend     = as_vector( [-1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*exp(-2.*pi*pi*mu*0.1), sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*mu*0.1), -0.25*(cos(2.*pi*x[0])+cos(2.*pi*x[1]))*exp(-4.*pi*pi*mu*time)])


    f           = as_vector( [2.*pi*pi*mu*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t , -2.*pi*pi*mu*sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*mu*time) ] )
    f           += as_vector( [0.25*2.*pi*sin(2.*pi*x[0])*exp(-4.*pi*pi*mu*time),  0.25*2.*pi*sin(2.*pi*x[1])*exp(-4.*pi*pi*mu*time)] )

    f           -= as_vector( [mu*2*gamma_t*pi*pi*cos(1.*pi*x[0])*sin(1.*pi*x[1]), mu*-2*gamma_t*pi*pi*sin(1.*pi*x[0])*cos(1.*pi*x[1])] )
    f           += as_vector( [cos(1.*pi*x[0])*sin(1.*pi*x[0])*(cos(1.*pi*x[1])*cos(1.*pi*x[1])-sin(1.*pi*x[1])*sin(1.*pi*x[1])),cos(1.*pi*x[1])*sin(1.*pi*x[1])*(cos(1.*pi*x[0])*cos(1.*pi*x[0])-sin(1.*pi*x[0])*sin(1.*pi*x[0]))] )

    # f          += nu*exact_u

    solution = spc.interpolate( [exact_u[0],exact_u[1],exact_u[2]], name="velocity")
    velocity = as_vector([solution[0],solution[1]])
    pressure = as_vector([solution[2]])

    vtk = grid.sequencedVTK("monolithic-Taylor-Green",
            pointvector={"velocity":velocity,"pressure":pressure})
    vtk()

    # get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
    # old_solution = solution.copy()
    old_sol = solution.copy(name="uOld")


    dirichlet = conditional(abs(x[0])<1e-8,1.,0.)

    du = 0.5 *(grad(u)+grad(u).T)

    # Formulation from Kroener Toulopoulos paper on P-NS
    # b = (pow(d + sqrt(inner(du, du)), (pnb-2))*inner(du, grad(v)) ) * dx

    abs_du = (inner(du, du))

    # GG = du[].two_norm
    b = (pow(d**(1) + (abs_du+1e-50)**0.5, (pnb-2)/1)*inner(du, grad(v)) ) * dx

    # b = (pow(d + inner(grad(u), grad(u)), (pnb-2)/2)*inner(grad(u), grad(v)) ) * dx


    tau = NamedConstant(spc,name="tau")
    a = (inner(u - as_vector([old_sol[0],old_sol[1]]), v)) * dx + tau * viscosity * b
         # tau * viscosity * inner(grad(u), grad(v))) * dx


    a -= tau * inner(p[0],div(v)) * dx
    a += tau * inner(grad(u), outer(v, u)) * dx
    a += tau * inner(div(u), q[0]) * dx
    a += 1e-3 * (inner((p[0] - old_sol[2]), q[0])) * dx
    a -= tau * inner(f,v) * dx

    model = create.model("elliptic", grid, a == 0,
            # DirichletBC(spc, [None,None,None], 2),       # bottom/top
            # DirichletBC(spc, [None,None,None], 3),       # bottom/top
            # DirichletBC(spc, [0.5,0,None], 4),             # left
            # DirichletBC(spc, [0,0,None], 1))             # right
            DirichletBC(spc,[exact_u[0],exact_u[1],exact_u[2]],1))
    model.tau = deltaT

    # setup structure for olver parameters
    solverParameter={"fem.solver.newton.absolutetol": 1e-13,
                     "fem.solver.newton.reductiontol": 1e-13,
                     "fem.solver.newton.tolerance": 1e-12,
                     "fem.solver.newton.verbose": "true",
                     "fem.solver.newton.linear.verbose": "false"}
    # create the solver using a standard fem scheme
    scheme = create.scheme("h1", model, spc, parameters=solverParameter)

    # mainOp.model.mu = 0.01
    # mainOp.model.nu = 0.0

    def solvestep1():
        old_sol.assign(solution)
        scheme.solve(target=solution)


    endTime =50
    timeStep = deltaT
    time = timeStep
    while time < endTime:
        print( "Time is:", time )
        solvestep1()
        vtk()
        l2error_fn = inner ( as_vector([solution[0],solution[1]])- as_vector([exact_uend[0],exact_uend[1]]), as_vector([solution[0],solution[1]]) - as_vector([exact_uend[0],exact_uend[1]]))
        l2error = sqrt( integrate(grid, l2error_fn, 5)[0] )
        print('|u_h - u| =', l2error)
        time += timeStep

compute()
