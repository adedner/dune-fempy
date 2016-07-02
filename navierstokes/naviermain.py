from __future__ import print_function
from mpi4py import MPI
import math
import ufl
import dune.models.femufl as duneuflmodel
import dune.fem as fem
import dune.fem.scheme as scheme

ufl2model = duneuflmodel.DuneUFLModel(2,3) # this is utility.init and sets the dim range
u = ufl2model.trialFunction()
v = ufl2model.testFunction()
x = ufl2model.spatialCoordinate()
a = (ufl.inner(u,v) + ufl.inner(ufl.grad(u),ufl.grad(v)))*ufl.dx(0)
bnd_u = ufl2model.coefficient("bnd_u",1)
g = [bnd_u,0,0]
ufl2model.generate(a,diric={1:g})

# initialise grid
grid2d = fem.leafGrid( "../data/hole2.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming" )
grid2d.hierarchicalGrid.globalRefine(6)
# grid2d = fem.leafGrid( (fem.reader.gmsh,"../data/karmanvortexstreet.msh"), "ALUCubeGrid", dimgrid=2 )
# grid2d.hierarchicalGrid.globalRefine(1)

viscosity = 0.003
timeStep = 0.005
endTime = 70
saveinterval = 0.1
problemNumber = 4   # this should be a model

# spaces
pressureSpace = fem.create.space( "Lagrange", grid2d, polorder = 1, dimrange = 1 )
velocitySpace = fem.create.space( "Lagrange", grid2d, polorder = 2, dimrange = 2 )
# velocitySpace = fem.create.space( "P1Bubble", grid2d, dimrange=2 )
# probem with missing dirichlet points in bubble space - need to update
# dirichletconstraints

# schemes
model = ufl2model.makeAndImport(grid2d,name="laplace").get()
stokesScheme = fem.create.scheme( "StokesScheme", ( velocitySpace, pressureSpace),model,\
               viscosity, problemNumber, "stokes", timeStep, storage = "Istl" )
burgersScheme = fem.create.scheme( "BurgersScheme", ( velocitySpace, pressureSpace),model,\
               viscosity, problemNumber, "burgers", timeStep, storage = "Istl" )

# set up solution initializating with data at t=0
velocity = velocitySpace.interpolate( lambda x: [0,0], name = "velocity", storage = "Istl" )
pressure = pressureSpace.interpolate( lambda x: [0], name = "pressure", storage = "Istl" )
solution = velocity, pressure

# vtk output of solution + voricity
def vorticityLocal(element,x):
    jac = velocity.localFunction(element).jacobian(x)
    return [ jac[1][0] - jac[0][1] ]
vorticity = grid2d.localGridFunction( "vorticity", vorticityLocal)
vtk = grid2d.writeVTK( "ns_", pointdata=solution+(vorticity,), number=0 )

# time loop
stokesScheme.next()     # to be removed (increases the time for forcing,bc)
burgersScheme.next()    # to be removed
time = timeStep
counter = 0
save_counter = 1
savestep = saveinterval
def bnd_u(x):
    if time < 0.05:
      return [ (1+x[0])*(1-x[1]) * time / 0.05 ]
    return [(1+x[0])*(1-x[1])]
bnd_uGlobal = grid2d.globalGridFunction("bnd_velocity", bnd_u)
model.setbnd_u( bnd_uGlobal )
while time < endTime:
    print( "Time is:", time )
    print( 'Solve step 1 - Stokes' )
    stokesScheme.solve( rhs = solution, target = solution, assemble = counter==0 )
    print( 'Solve step 2 - Burgers' )
    burgersScheme.solve( rhs = solution, target = solution, assemble = counter==0 )
    print( 'Solve step 3 - Stokes' )
    stokesScheme.solve( rhs = solution, target = solution, assemble = False )
    time += timeStep
    stokesScheme.next()   # to be removed
    burgersScheme.next()  # to be removed
    if time > savestep:
        savestep += saveinterval
        vtk.write( "ns_", save_counter )
        save_counter += 1
    counter += 1
