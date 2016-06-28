from mpi4py import MPI
import dune.fem as fem
import dune.fem.scheme as scheme

# initialise grid
grid2d = fem.leafGrid( "../data/hole2_larger.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming" )
#grid2d = grid.leafGrid( "../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2 )

grid2d.hierarchicalGrid.globalRefine(1)

timeStep = 0.001
endTime = 10
problemNumber = 4
velocitySpace = fem.create.space( "Lagrange", grid2d, polorder = 2, dimrange = 2 )
pressureSpace = fem.create.space( "Lagrange", grid2d, polorder = 1, dimrange = 1 )
stokesScheme = fem.create.scheme( "StokesScheme", ( velocitySpace, pressureSpace),\
               problemNumber, "stokes", timeStep, storage = "Istl" )
burgersScheme = fem.create.scheme( "BurgersScheme", ( velocitySpace, pressureSpace),\
               problemNumber, "stokes", timeStep, storage = "Istl" )
velocity = velocitySpace.interpolate( lambda x: [0,0], name = "velocity", storage = "Istl" )
pressure = pressureSpace.interpolate( lambda x: [0], name = "pressure", storage = "Istl" )
solution = velocity, pressure

stokesScheme.initialize( solution )

vtk = grid2d.writeVTK( "ns_", pointdata = solution, number = 0 )

def solve_method( timeStep, endTime ):
    stokesScheme.next()
    burgersScheme.next()
    time = timeStep
    counter = 0
    while time < endTime:
        print( "Time is:", time )
        print( 'Solve step 1 - Stokes' )
        stokesScheme.solve( rhs = solution, target = solution, assemble = counter==0 )
        print( 'Solve step 2 - Burgers' )
        burgersScheme.solve( rhs = solution, target = solution, assemble = counter==0 )
        print( 'Solve step 3 - Stokes' )
        stokesScheme.solve( rhs = solution, target = solution, assemble = False )
        time += timeStep
        stokesScheme.next()
        burgersScheme.next()
        counter += 1
        vtk.write( "ns_", counter )

solve_method( timeStep, endTime )

print( 'Finished' )
