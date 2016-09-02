from __future__ import print_function
from mpi4py import MPI
import math
from ufl import *
from dune.ufl import Space as UFLSpace
from dune.models.elliptic import importModel
import dune.fem as fem
from dune.femmpi import parameter

parameter.append("../data/parameter-navier")

# initialise grid
grid = fem.leafGrid( "../data/hole2.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming" )
# grid = fem.leafGrid( (fem.reader.gmsh,"../data/karmanvortexstreet.msh"), "ALUCubeGrid", dimgrid=2 )
grid.hierarchicalGrid.globalRefine(6)

viscosity = 0.03
timeStep = 0.05
endTime = 70

# boundary condition
saveinterval = 0.1
def inflow_u(x):
    ux = 0
    if x[0]<-1+1e-8:
        ux = min(1.0,(((x[1]+1.)*(1.-x[1])*time)/(10.*timeStep)))
    return [ux,0,0]

# model
uflSpace    = UFLSpace((2,2), 3)
u           = TrialFunction(uflSpace)
v           = TestFunction(uflSpace)
x           = SpatialCoordinate(uflSpace.cell())
bnd_u       = Coefficient(uflSpace)

a = inner(grad(u),grad(v)) * dx(0)
model = importModel(grid, a == 0, dirichlet={1:bnd_u}, tempVars=False).get()
model.setCoefficient(bnd_u, grid.globalGridFunction("bnd", inflow_u))

# spaces
pressureSpace = fem.create.space( "Lagrange", grid, polorder = 1, dimrange = 1 )
velocitySpace = fem.create.space( "Lagrange", grid, polorder = 2, dimrange = grid.dimWorld )
# velocitySpace = fem.create.space( "P1Bubble", grid, # dimrange=grid.dimWorld )
# problem with missing dirichlet points in bubble space - need to update
# dirichletconstraints

# schemes
stokesScheme = fem.create.scheme( "StokesScheme", ( velocitySpace, pressureSpace), model, "stokes",\
               viscosity, timeStep, storage = "Istl" )
burgersScheme = fem.create.scheme( "BurgersScheme", ( velocitySpace, pressureSpace), model, "burgers",\
                viscosity, timeStep, storage = "Istl" )

# set up solution initializating with data at t=0
velocity = velocitySpace.interpolate( lambda x: [0,0], name = "velocity", storage = "Istl" )
pressure = pressureSpace.interpolate( lambda x: [0], name = "pressure", storage = "Istl" )
solution = velocity, pressure

# vtk output of solution + voricity
def vorticityLocal(element,x):
    jac = velocity.localFunction(element).jacobian(x)
    return [ jac[1][0] - jac[0][1] ]
vorticity = grid.localGridFunction( "vorticity", vorticityLocal)
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
