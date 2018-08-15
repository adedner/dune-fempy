from __future__ import print_function
from mpi4py import MPI
import math

from ufl import *
from dune.ufl import Space as UFLSpace

import dune.create as create
import dune.fem

dune.fem.parameter.append("../data/parameter-navier")

# initialise grid
#grid = create.grid( "ALUSimplex", "../data/hole2.dgf", dimgrid=2, refinement="conforming" )
grid = create.grid( "ALUSimplex", "../data/hole2.dgf", dimgrid=2 )
grid.hierarchicalGrid.globalRefine(6)

viscosity = 0.03
timeStep = 0.05
endTime = 70

# boundary condition
saveinterval = 0.1
def inflow_u(x):
    ux = 0
    if x[0]<-1+1e-8:
        ux = (x[1]+1.)*(1.-x[1]) * min(time/(10.*timeStep),1)
    return [ux,0,0]

# model
uflSpace    = UFLSpace((2,2), 3)
u           = TrialFunction(uflSpace)
v           = TestFunction(uflSpace)
x           = SpatialCoordinate(uflSpace.cell())
bnd_u       = Coefficient(uflSpace)

a = inner(grad(u),grad(v)) * dx(0)
bnd = create.function("global", grid, "bnd", 3, inflow_u)
model = create.model("elliptic", grid, a == 0)

# spaces
pressureSpace = create.space( "Lagrange", grid, polorder = 1, dimrange = 1, storage="istl" )
# velocitySpace = fem.space.create( "Lagrange", grid, polorder = 2, dimrange = grid.dimWorld )
velocitySpace = create.space( "P1Bubble", grid, dimrange=grid.dimWorld, storage="istl" )
# problem with missing dirichlet points in bubble space - need to update
# dirichletconstraints
# set up initial conditions
# solution = spc.interpolate(lambda x: [math.atan((10.0 * x[0] * (1-x[0]) * x[1] * (1-x[1]))**2)], name="u")
# vtk = grid.sequencedVTK("heat", pointdata=[solution])
# vtk()
# schemes
# stokesScheme = create.scheme( "stokes", ( velocitySpace, pressureSpace ), model, "stokes", viscosity, timeStep )
# burgersScheme = create.scheme( "burgers", ( velocitySpace, pressureSpace ), model, "burgers",\
#                 viscosity, timeStep )

# set up solution initializating with data at t=0
velocity = velocitySpace.interpolate( lambda x: [ x[1] * ( 1.0 - x[ 1 ] ),0], name = "velocity", storage = "Istl" )
pressure = pressureSpace.interpolate( lambda x: [0], name = "pressure", storage = "Istl" )
solution = velocity, pressure


# def burgerScheme(oldsolution,timestep):
#     u_n = oldsolution
u_n = velocity.copy();
burg = (inner(u - u_n, v) + deltaT * inner(grad(u), outer(v, u))) * dx
# a = (inner(u - u_n, v) + inner(grad(u), outer(v, u_n))) * dx


model = create.model("integrands", grid, burg == 0)

# setup structure for olver parameters
solverParameter={"fem.solver.newton.linabstol": 1e-11,
                 "fem.solver.newton.linreduction": 1e-11,
                 "fem.solver.newton.tolerance": 1e-10,
                 "fem.solver.newton.verbose": "true",
                 "fem.solver.newton.linear.verbose": "false"}
viscosity = 0.03
timeStep = 0.05
# create the solver using a standard fem scheme
# scheme = create.scheme("h1", spc, model, parameters=solverParameter)
# scheme = create.scheme("h1galerkin", spcP, model, parameters=solverParameter)
#scheme = create.scheme("dggalerkin", spcP, model, 15*theta*deltaT, parameters=solverParameter)
# burgersScheme = create.scheme( "burgers", ( spcU, spcP ), model, "burgers",\
#            viscosity, deltaT )
# scheme = create.scheme("galerkin", model, spcP, parameters=solverParameter)
scheme = create.scheme("galerkin", burg == 0, spcU, solver="gmres", parameters = solverParameter)

old_solution.assign(velocity)
scheme.solve(target=velocity)

# vtk = grid.sequencedVTK("heat", pointdata=[solution])
# vtk()
# vtk output of solution + voricity
# def vorticityLocal(element,x):
#     jac = velocity.localFunction(element).jacobian(x)
#     return [ jac[1][0] - jac[0][1] ]
# vorticity = grid.localGridFunction( "vorticity", 2)
def plot(count=None):
    grid.writeVTK("Stokes",
            pointdata={"pressure":pressure},
            pointvector={"velocity":velocity},
            number=count
    )

# vtk = grid.writeVTK( "ns_", pointdata=solution, number=0 )

plot()
# time loop
time = timeStep
counter = 0
savestep = saveinterval
while time < endTime:
    print( "Time is:", time )
    print( 'Solve step 1 - Stokes' )
    burgerScheme(velocity,timeStep)
    # stokesScheme.solve( rhs = solution, target = solution, assemble = counter==0 )
    print( 'Solve step 2 - Burgers' )
    # burgersScheme.solve( rhs = solution, target = solution, assemble = counter==0 )
    print( 'Solve step 3 - Stokes' )
    # stokesScheme.solve( rhs = solution, target = solution, assemble = False )
    time += timeStep
    if time > savestep:
        savestep += saveinterval
        vtk.write( "ns_", counter + 1 )
    counter += 1
