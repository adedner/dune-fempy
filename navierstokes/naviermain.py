from __future__ import print_function
from mpi4py import MPI
import math

import dune
dune.initialize()

from ufl import *
from dune.ufl import Space as UFLSpace
# from dune.models.elliptic import importModel
import dune.fem as fem
from dune.femmpi import parameter
import dune.grid
import dune.alugrid
import dune.fem.space
import dune.fem.scheme

parameter.append("../data/parameter-navier")

# initialise grid
grid = dune.grid.create( "ALUSimplex", "../data/hole2.dgf", dimgrid=2, refinement="conforming" )
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
model = dune.fem.create.ellipticModel(grid, a == 0, dirichlet={1:bnd_u})\
        ( coefficients={ bnd_u: grid.globalGridFunction("bnd", 3, inflow_u) } )

# spaces
pressureSpace = fem.space.create( "Lagrange", grid, polorder = 1, dimrange = 1 )
# velocitySpace = fem.space.create( "Lagrange", grid, polorder = 2, dimrange = grid.dimWorld )
velocitySpace = fem.space.create( "P1Bubble", grid, dimrange=grid.dimWorld )
# problem with missing dirichlet points in bubble space - need to update
# dirichletconstraints

# schemes
stokesScheme = fem.scheme.create( "Stokes", ( velocitySpace, pressureSpace ), model, "stokes",\
               viscosity, timeStep, storage = "Istl" )
burgersScheme = fem.scheme.create( "Burgers", ( velocitySpace, pressureSpace ), model, "burgers",\
                viscosity, timeStep, storage = "Istl" )

# set up solution initializating with data at t=0
velocity = velocitySpace.interpolate( lambda x: [0,0], name = "velocity", storage = "Istl" )
pressure = pressureSpace.interpolate( lambda x: [0], name = "pressure", storage = "Istl" )
solution = velocity, pressure

# vtk output of solution + voricity
def vorticityLocal(element,x):
    jac = velocity.localFunction(element).jacobian(x)
    return [ jac[1][0] - jac[0][1] ]
vorticity = grid.localGridFunction( "vorticity", 2, vorticityLocal)
vtk = grid.writeVTK( "ns_", pointdata=solution+(vorticity,), number=0 )

# time loop
time = timeStep
counter = 0
savestep = saveinterval
while time < endTime:
    print( "Time is:", time )
    print( 'Solve step 1 - Stokes' )
    stokesScheme.solve( rhs = solution, target = solution, assemble = counter==0 )
    print( 'Solve step 2 - Burgers' )
    burgersScheme.solve( rhs = solution, target = solution, assemble = counter==0 )
    print( 'Solve step 3 - Stokes' )
    stokesScheme.solve( rhs = solution, target = solution, assemble = False )
    time += timeStep
    if time > savestep:
        savestep += saveinterval
        vtk.write( "ns_", counter + 1 )
    counter += 1
