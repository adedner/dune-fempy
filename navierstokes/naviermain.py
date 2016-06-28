from mpi4py import MPI
import dune.fem as fem
import dune.fem.scheme as scheme

# initialise grid
grid2d = fem.leafGrid( (fem.reader.gmsh,"../data/karmanvortexstreet.msh"), "ALUCubeGrid", dimgrid=2 )
# grid2d.hierarchicalGrid.globalRefine(1)

viscosity = 0.01
timeStep = 0.005
endTime = 10
saveinterval = 0.1
problemNumber = 4   # this should be a model

# spaces
pressureSpace = fem.create.space( "Lagrange", grid2d, polorder = 1, dimrange = 1 )
velocitySpace = fem.create.space( "Lagrange", grid2d, polorder = 2, dimrange = 2 )
# velocitySpace = fem.create.space( "P1Bubble", grid2d, dimrange=2 )
# probem with missing dirichlet points in bubble space - need to update
# dirichletconstraints

# schemes
stokesScheme = fem.create.scheme( "StokesScheme", ( velocitySpace, pressureSpace),\
               viscosity, problemNumber, "stokes", timeStep, storage = "Istl" )
burgersScheme = fem.create.scheme( "BurgersScheme", ( velocitySpace, pressureSpace),\
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
