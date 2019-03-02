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

# Crank Nicholson
# theta = 1.
deltaT = 1e-2 # 0.001
viscosity = 1.0e-2
d = 0.001
pnb = 2.5
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
    grid = create.view("adaptive", grid="ALUCube", constructor="../../data/channelexpansion.dgf", dimgrid=2)
    grid.hierarchicalGrid.globalRefine(6)
    # set up a lagrange scalar space with polynomial order 2 over that grid
    spc = create.space("lagrange", grid, dimrange=3, order=2, storage="istl")

    # set up initial conditions

    solution = spc.interpolate( [0,0,0], name="velocity")
    velocity = as_vector([solution[0],solution[1]])
    pressure = as_vector([solution[2]])

    vtk = grid.sequencedVTK("p2.5-coupledNavstokes",
            pointvector={"velocity":velocity,"pressure":pressure})
    vtk()

    # get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
    # old_solution = solution.copy()
    old_sol = solution.copy(name="uOld")

    # now define the actual pde to solve:
    #            u - u_n deltaT laplace( theta u + (1-theta) u_n ) = 0

    trial = TrialFunction(spc)
    test  = TestFunction(spc)
    u = as_vector([trial[0],trial[1]])
    v = as_vector([test[0],test[1]])
    p = as_vector([trial[2]])
    q = as_vector([test[2]])
    x = SpatialCoordinate(spc)



    def mark(element):
        solutionLocal = old_sol.localFunction(element)
        grad = solutionLocal.evaluate(element.geometry.referenceElement.center)
        if grad[0] > 0.5:
            return Marker.refine if element.level < maxLevel else Marker.keep
        else:
             return Marker.keep


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
    a += 400 * (inner((p[0] - old_sol[2]), q[0])) * dx

    model = create.model("elliptic", grid, a == 0,
            DirichletBC(spc, [None,None,None], 2),       # bottom/top
            DirichletBC(spc, [None,None,None], 3),       # bottom/top
            DirichletBC(spc, [0.5,0,None], 4),             # left
            DirichletBC(spc, [0,0,None], 1))             # right

    model.tau = deltaT

    # setup structure for olver parameters
    solverParameter={"fem.solver.newton.absolutetol": 1e-13,
                     "fem.solver.newton.reductiontol": 1e-13,
                     "fem.solver.newton.tolerance": 1e-12,
                     "fem.solver.newton.verbose": "true",
                     "fem.solver.newton.linear.verbose": "false"}
    # create the solver using a standard fem scheme
    scheme = create.scheme("h1", model, spc, parameters=solverParameter)

    endTime =50
    timeStep = deltaT
    time = timeStep
    while time < endTime:
        print( "Time is:", time )
        old_sol.assign(solution)

        # grid.hierarchicalGrid.mark(mark)
        # fem.adapt(grid.hierarchicalGrid,[solution])
        # fem.loadBalance(grid.hierarchicalGrid,[solution])

        scheme.solve(target=solution)

        vtk()
        time += timeStep

compute()
