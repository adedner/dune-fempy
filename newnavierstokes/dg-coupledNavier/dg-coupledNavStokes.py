from __future__ import print_function

import math
from ufl import *

from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, dot, div, grad, dx, as_vector, transpose, Identity

from dune.grid import cartesianDomain, Marker

from dune.ufl import NamedConstant

import dune.create as create

from dune.fem.function import levelFunction

from dune.fem import parameter, adapt

from dune.ufl import DirichletBC, Space

import dune.fem as fem

# Crank Nicholson
theta = 1
deltaT = 1e-2# 0.001
viscosity = 1.0e-2
mup = 1e1
mu  = 1e3
d = 0.001
pnb = 2.
maxLevel = 5
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

    grid.hierarchicalGrid.globalRefine(4)
    # set up a lagrange scalar space with polynomial order 2 over that grid
    spc = create.space("dgonbhp", grid, dimrange=3, order=3, storage="istl")
    # spcP = create.space("lagrange", grid, dimrange=1, order=2)
    spcP = create.space("dgonb", grid, dimrange=1, order=1, storage="istl")


    fvspc = create.space("finitevolume", grid, dimrange=3, storage="istl")
    estimate = fvspc.interpolate([0,0,0], name="estimate")
    # set up initial conditions

    solution = spc.interpolate( [0,0,0], name="velocity")
    velocity = as_vector([solution[0],solution[1]])
    pressure = as_vector([solution[2]])


    # velocity = spc.interpolate( [0,0], name="velocity")
    # pressure = spcP.interpolate( lambda x: [0], name = "pressure", storage = "Istl" )

    # solution = velocity, pressure
    #
    # vtk = grid.sequencedVTK("coupledNavstokes", pointdata={"velocity":velocity})
    # vtk()

    vtk = grid.sequencedVTK("dg-coupledNavstokes",
            pointvector={"velocity":velocity,"pressure":pressure})
    vtk()

    # get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
    old_sol = solution.copy()
    # old_velo = velocity.copy()

    # old_solution.name = "uOld"
    # old_velo.name = "uOld"

    # now define the actual pde to solve:
    #            u - u_n deltaT laplace( theta u + (1-theta) u_n ) = 0
    # uflSpace = Space((grid.dimGrid, grid.dimWorld), 2)
    # uflSpaceP = Space((grid.dimGrid, grid.dimWorld), 2)

    trial = TrialFunction(spc)
    test  = TestFunction(spc)
    u = as_vector([trial[0],trial[1]])
    v = as_vector([test[0],test[1]])
    p = as_vector([trial[2]])
    q = as_vector([test[2]])
    x = SpatialCoordinate(spc)

    # u = TrialFunction(uflSpace)
    # v = TestFunction(uflSpace)
    # p = TrialFunction(uflSpaceP)
    # q = TestFunction(uflSpaceP)
    # x = SpatialCoordinate(uflSpace.cell())
    # n = FacetNormal(spc.cell())
    # n, h = FacetNormal(spc.cell()), MinFacetEdgeLength(spc.cell())
    hT = MaxCellEdgeLength(spc.cell())
    he = MaxFacetEdgeLength(spc.cell())('+')
    n = FacetNormal(spc.cell())
    # mu = 1000/ avg(h)
    # mu = 7.5 * 16

    estimator_ufl = hT**2 * (div(grad(u[0])))**2 * v[0] * dx + he * inner(jump(grad(u[0])), n('+'))**2 * avg(v[0]) * dS
    estimator_model = create.model("integrands", grid, estimator_ufl == 0)
    estimator = create.operator("galerkin", estimator_model, spc, fvspc)

    tolerance = 0.001
    gridSize = grid.size(0)
    def mark(element):
        solutionLocal = solution.localFunction(element)
        grad = solutionLocal.jacobian(element.geometry.referenceElement.center)
        if grad[0].infinity_norm > 1.2:
            return Marker.refine if element.level < maxLevel else Marker.keep
        else:
            return Marker.coarsen
    # def mark(element):
    #     estLocal = estimate(element, element.geometry.referenceElement.center)
    #     return Marker.refine if estLocal[0] > tolerance / gridSize else Marker.keep

    phi = conditional(x[0] < 1e-2, 1, 0)


    du = 0.5 *(grad(u)+grad(u).T)

    # Formulation from Kroener Toulopoulos paper on P-NS
    # b = (pow(d + sqrt(inner(du, du)), (pnb-2))*inner(du, grad(v)) ) * dx

    abs_du = (inner(du, du))

    # GG = du[].two_norm
    b = (pow(d**(1) + (abs_du+1e-50)**0.5, (pnb-2)/1)*inner(du, grad(v)) ) * dx

    # b = (pow(d + inner(grad(u), grad(u)), (pnb-2)/2)*inner(grad(u), grad(v)) ) * dx


    tau = NamedConstant(spc,name="tau")
         # tau * viscosity * inner(grad(u), grad(v))) * dx

    u_d = as_vector( [phi, 0] )


    dirichlet = conditional(x[0]<1e-2,1.,0. ) + conditional(x[0]>29.99 ,1.,0. )

    # u_n = old_velo
    # tau = Constant(triangle)
    # a = (inner(u - u_n, v) + tau * viscosity * inner(grad(theta*u + (1-theta)*u_n), grad(v))) * dx
    # a += -1*tau*inner(p[0],div(v)) * dx
    # a += tau* inner(grad(u), outer(v, u)) * dx
    # a += div(u)*q[0] * dx


    # tau = NamedConstant(spc,"tau")

    # diffusion
    a = (inner(u - as_vector([old_sol[0],old_sol[1]]), v)) * dx + tau * viscosity * b

    # a = (inner(u - as_vector([old_sol[0],old_sol[1]]), v) +
    #          tau * viscosity * 0.5 * inner(grad(u)+grad(u).T, grad(v))) * dx
    # a  = tau * viscosity * inner(grad(u), grad(v)) * dx
    a -= tau * viscosity * (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(du), outer(jump(v), n('+')))) * dS
    a += tau * mu *  viscosity * inner(jump(u), jump(v)) * dS
    # a -= tau * (inner(outer(u_d, n), grad(v)) + inner(grad(u_d), outer(v, n))) * ds
    a +=tau *  viscosity * (inner(outer(u_d, n), grad(v)) ) * dirichlet *  ds
    a += tau * viscosity *( mu /  (500))* inner(u_d, v) * dirichlet * ds

    # advection
    a += tau* 0.5 * inner(grad(u)+grad(u).T, outer(v, u)) * dx
    a += 0.5 * tau* inner(v, outer(div(u), u)) * dx
    a -= tau * (inner(outer(jump(u), n('+')), outer(avg(u),avg(v)) )) * dS
    # a -= 0.5 * tau * (inner(outer(jump(u), n('+')), avg(u*v) )) * dS

    # pressure part
    a += -1*tau*inner(p[0],div(v)) * dx
    a += tau * (inner(outer(avg(p[0]), n('+')), jump(v))) * dS
    # # a += tau * inner(p_d, v) * dirichlet * ds

    # # incompressibilty (pressure part)
    # # a -= inner(div(u),q[0]) * dx
    a += 1 * tau * (inner(outer(avg(q[0]), n('+')), jump(u))) * dS
    a += 1 * tau * (inner(outer(q[0], n), u_d)) * dirichlet * ds
    # # a += tau * (1/viscosity) * mup * inner(jump(u), outer(jump(q[0]), n('+'))) * dS
    # # a += tau * (1/viscosity) * mup * inner(u_d, outer(q[0], n)) *dirichlet *  ds
    #
    # # complete coupled scheme
    # # a = inner(u - old_velo, v) * dx + 1*theta*a + 1*(1-theta)*replace(a,{u:old_velo})
    #
    #
    # # a -= (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
    # # a += mu * inner(jump(u), jump(v)) * dS
    # # a -= (inner(outer(u, n), grad(v)) + inner(grad(u), outer(v, n))) * ds
    # # a += mu * inner(u, v) * ds
    # # now generate the model code and compile
    # #model = create.model("elliptic", grid, a == 0, coefficients={u_n:old_solution})
    # #model.setConstant(tau,[deltaT])

    a += -1 *tau * inner(div(u), q[0]) * dx

    # a += 1*1e-1*inner(p[0], q[0]) * dx

    a += 400 * (inner((p[0] - old_sol[2]), q[0])) * dx

    model = create.model("integrands", grid, a == 0)
    # model = create.model("elliptic", grid, a == 0,
    #         DirichletBC(uflSpace, [None,None], 2),       # bottom
    #         DirichletBC(uflSpace, [None,None], 3),      # top
    #         DirichletBC(uflSpace, [1,0], 4),          # left
    #         DirichletBC(uflSpace, [0,0], 1))                # right

    model.setConstant(tau, deltaT)

    # setup structure for olver parameters
    solverParameter={"fem.solver.newton.linabstol": 1e-7,
                     "fem.solver.newton.linreduction": 1e-7,
                     "fem.solver.newton.tolerance": 1e-6,
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
    # steps = 5000
    # for n in range(1,steps+1):
    #     # old_solution.assign(solution)
    #     old_velo.assign(velocity)
    #
    #     scheme.solve(target=velocity)
    #     vtk()
    endTime =50
    timeStep = deltaT
    time = timeStep
    while time < endTime:
        print( "Time is:", time )
        old_sol.assign(solution)
        scheme.solve(target=solution)

        # grid.hierarchicalGrid.mark(mark)
        # fem.adapt(grid.hierarchicalGrid,[solution])
        # fem.loadBalance(grid.hierarchicalGrid,[solution])

        vtk()
        time += timeStep

compute()
