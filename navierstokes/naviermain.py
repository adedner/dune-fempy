from mpi4py import MPI
import dune.fem as fem
import dune.fem.scheme as scheme

# initialise grid
grid2d = fem.leafGrid( (fem.reader.gmsh,"../data/karmanvortexstreet.msh"), "ALUSimplexGrid", dimgrid=2 )

grid2d.hierarchicalGrid.globalRefine(1)

timeStep = 0.005
endTime = 10
saveinterval = 0.1
problemNumber = 4
velocitySpace = fem.create.space( "Lagrange", grid2d, polorder=2, dimrange=2 )
# velocitySpace = fem.create.space( "P1Bubble", grid2d, dimrange=2 )
# probem with missing dirichlet points in bubble space - need to update
# dirichletconstraints
pressureSpace = fem.create.space( "Lagrange", grid2d, polorder=1, dimrange=1 )
stokesScheme  = fem.create.scheme( "StokesScheme", ( velocitySpace, pressureSpace),\
                           problemNumber,"stokes", timeStep )
bs = scheme.get( "BurgersScheme", ( velocitySpace, pressureSpace ) )
burgersScheme = bs.Scheme( ( velocitySpace, pressureSpace ), problemNumber, timeStep )

velocitySpace = 0
pressureSpace = 0

stokesScheme.initialize()

def vorticityLocal(element,x):
    jac = stokesScheme.solution()[0].localFunction(element).jacobian(x)
    return [ jac[1][0] - jac[0][1] ]
vorticity = grid2d.localGridFunction( "vorticity", vorticityLocal)

vtk = grid2d.writeVTK( "ns_", pointdata=stokesScheme.solution()+(vorticity,), number=0 )

def solve_method( timeStep, endTime ):
    stokesScheme.next()
    burgersScheme.next()
    time = timeStep
    counter = 0
    save_counter = 1
    savestep = saveinterval
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
        if time > savestep:
            savestep += saveinterval
            vtk.write( "ns_", save_counter )
            save_counter += 1
        counter += 1

solve_method( timeStep, endTime )

print( 'Finished' )
