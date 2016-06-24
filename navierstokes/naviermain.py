from mpi4py import MPI

import dune.fem.function as gf
import dune.fem.grid as grid
import dune.fem.scheme as scheme
import dune.fem.space as space

# initialise grid
grid2d = grid.leafGrid( "../data/hole2_larger.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming" )
#grid2d = grid.leafGrid( "../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2 )

grid2d.hierarchicalGrid.globalRefine(3)

timeStep = 0.001
endTime = 10
problemNumber = 4
# initialise Stokes scheme
space = space.create( "Lagrange", grid2d ) # not actually used
ss = scheme.get( "StokesScheme", space, grid2d, 1 ) # ideally this should be ss = scheme.get( "StokesScheme" )
stokesScheme = ss.StokesScheme( grid2d, problemNumber, timeStep )
# initialise Burgers scheme
bs = scheme.get( "BurgersScheme", space, grid2d, 1 )
burgersScheme = bs.BurgersScheme( grid2d, problemNumber, timeStep )
stokesScheme.initialize()

vtk = grid2d.vtkWriter()
stokesScheme.solution()[0].addToVTKWriter( vtk, vtk.PointData )
stokesScheme.solution()[1].addToVTKWriter( vtk, vtk.PointData )

def solve_method( timeStep, endTime ):
    counter = 0
    vtk.write( "ns_0000" )
    stokesScheme.next()
    burgersScheme.next()
    time = timeStep
    while time < endTime:
        print( time, " burgers=", burgersScheme.time(), "  stokes=", stokesScheme.time())
        # first step (solve Stokes for velocity and pressure)
        print( 'Solve step 1 - Stokes' )
        stokesScheme.preparestep1()
        stokesScheme.solve( counter == 0 ) # assemble only on first iteration
        # stokesScheme.solve( rhs=burgerScheme.solution(), target=stokesScheme.solution(), counter==0 )

        # second step (solve Burgers for velocity)
        burgersScheme.update( stokesScheme.solution() )
        print( 'Solve step 2 - Burgers' )
        burgersScheme.prepare()
        burgersScheme.solve( counter == 0 )
        # burgersScheme.solve( rhs=stokesScheme.solution(), target=burgersScheme.solution(), counter==0 )

        # third step (solve Stokes for velocity and pressure)
        stokesScheme.update( burgersScheme.solution() )
        print( 'Solve step 3 - Stokes' )
        stokesScheme.preparestep3()
        stokesScheme.solve( False )
        # stokesScheme.solve( rhs=burgerScheme.solution(), target=stokesScheme.solution(), False )

        time += timeStep
        stokesScheme.next()
        burgersScheme.next()
        counter += 1
        if abs(time%0.1) < 0.001:
            outName = str( counter ).zfill( 4 )
            vtk.write( "ns_"+outName )

solve_method( timeStep, endTime )

print( 'Finished' )
