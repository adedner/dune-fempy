from mpi4py import MPI
import dune.fem as fem
import dune.fem.scheme as scheme

# initialise grid
grid2d = fem.leafGrid( "../data/hole2_larger.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming" )
#grid2d = grid.leafGrid( "../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2 )

grid2d.hierarchicalGrid.globalRefine(5)

timeStep = 0.001
endTime = 10
problemNumber = 4
velocitySpace = fem.create.space( "Lagrange", grid2d, polorder=2, dimrange=2 )
pressureSpace = fem.create.space( "Lagrange", grid2d, polorder=1, dimrange=1 )
# ss = scheme.get( "StokesScheme", ( velocitySpace, pressureSpace) )
# stokesScheme = ss.Scheme( ( velocitySpace, pressureSpace ), problemNumber, timeStep )
stokesScheme = fem.create.scheme( "StokesScheme", ( velocitySpace, pressureSpace),\
               problemNumber,"stokes", timeStep )
bs = scheme.get( "BurgersScheme", ( velocitySpace, pressureSpace ) )
                   # velocitySpace=velocitySpace, pressureSpace=pressureSpace
burgersScheme = bs.Scheme( ( velocitySpace, pressureSpace ), problemNumber, timeStep )

velocitySpace = 0
pressureSpace = 0

stokesScheme.initialize()

vtk = grid2d.writeVTK( "ns_", pointdata=stokesScheme.solution(), number=0 )

def solve_method( timeStep, endTime ):
    stokesScheme.next()
    burgersScheme.next()
    time = timeStep
    counter = 0
    while time < endTime:
        print( "Time is:", time )
        print( 'Solve step 1 - Stokes' )
        stokesScheme.solve( rhs = burgersScheme.solution(), target = stokesScheme.solution(), assemble = counter==0 )
        print( 'Solve step 2 - Burgers' )
        burgersScheme.solve( rhs = stokesScheme.solution(), target = burgersScheme.solution(), assemble = counter==0 )
        print( 'Solve step 3 - Stokes' )
        stokesScheme.solve( rhs = burgersScheme.solution(), target = stokesScheme.solution(), assemble = False )
        time += timeStep
        stokesScheme.next()
        burgersScheme.next()
        counter += 1
        vtk.write( "ns_", counter )

solve_method( timeStep, endTime )

print( 'Finished' )
