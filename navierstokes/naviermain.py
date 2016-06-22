from mpi4py import MPI

import dune.fem.function as gf
import dune.fem.grid as grid
import dune.fem.scheme as scheme
import dune.fem.space as space

# initialise grid
grid2d = grid.leafGrid( "../data/hole2_larger.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming" )
#grid2d = grid.leafGrid( "../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2 )

#grid2d.hierarchicalGrid.globalRefine(2)

timeStep = 0.01
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
stokesScheme.pressure().addToVTKWriter( vtk, vtk.PointData )
stokesScheme.velocity().addToVTKWriter( vtk, vtk.PointData )

def solve_method( timeStep, endTime ):
    time = 0
    counter = 0
    outName = str( counter ).zfill( 4 )
    vtk.write( "ns_"+outName )
    while time < endTime:
        # first step (solve Stokes for velocity and pressure)
        print( 'Solve step 1 - Stokes' )
        stokesScheme.preparestep1()
        stokesScheme.solve( time == 0 ) # assemble only on first iteration

        # second step (solve Burgers for velocity)
        burgersScheme.updatepressure( stokesScheme.pressure() )
        burgersScheme.updatevelocity( stokesScheme.velocity() )
        print( 'Solve step 2 - Burgers' )
        burgersScheme.prepare()
        burgersScheme.solve( time == 0 )

        # third step (solve Stokes for velocity and pressure)
        stokesScheme.updatevelocity( burgersScheme.solution() )
        print( 'Solve step 3 - Stokes' )
        stokesScheme.preparestep3()
        stokesScheme.solve( False )

        time = time + timeStep
        stokesScheme.next()
        burgersScheme.next()
        counter = counter + 1
        outName = str( counter ).zfill( 4 )
        vtk.write( "ns_"+outName )

solve_method( timeStep, endTime )

print( 'Finished' )
