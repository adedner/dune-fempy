from mpi4py import MPI

import dune.fem.function as gf
import dune.fem.grid as grid
import dune.fem.scheme as scheme
import dune.fem.space as space

# initialise grid
grid2d = grid.leafGrid( "../data/hole2_larger.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming" )
#grid2d = grid.leafGrid( "../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2 )

grid2d.hierarchicalGrid.globalRefine(2)

timeStep = 0.001
endTime = 10
problemNumber = 4
velocitySpace = space.create( "Lagrange", grid2d, polorder=2, dimrange=2 )
pressureSpace = space.create( "Lagrange", grid2d, polorder=1 )
ss = scheme.get( "StokesScheme", ( velocitySpace, pressureSpace) )
stokesScheme = ss.Scheme( ( velocitySpace, pressureSpace ), problemNumber, timeStep )
bs = scheme.get( "BurgersScheme", ( velocitySpace, pressureSpace ) )
                   # velocitySpace=velocitySpace, pressureSpace=pressureSpace
burgersScheme = bs.Scheme( ( velocitySpace, pressureSpace ), problemNumber, timeStep )
stokesScheme.initialize()

vtk = grid2d.vtkWriter()
stokesScheme.solution()[0].addToVTKWriter( vtk, vtk.PointData )
stokesScheme.solution()[1].addToVTKWriter( vtk, vtk.PointData )

def solve_method( timeStep, endTime ):
    stokesScheme.next()
    burgersScheme.next()
    time = timeStep
    counter = 0
    vtk.write( "ns_0000" )
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
        # if abs(time%0.01) < 0.001:
        outName = str( counter ).zfill( 4 )
        vtk.write( "ns_"+outName )

solve_method( timeStep, endTime )

print( 'Finished' )
