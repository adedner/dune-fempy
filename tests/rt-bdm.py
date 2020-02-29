# NOTE: there is some issue with failing convergence when using solver=cg -
# it should work...
from dune.grid import structuredGrid, cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, dot, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import Constant, DirichletBC
from dune.fem.space import raviartThomas  as rtSpace
from dune.fem.space import bdm  as bdm

from dune.fem import parameter
parameter.append({"fem.verboserank": 0})

order = 1
# Note: structuredGrid fails in precon step!
# grid = structuredGrid([0,0],[3,1],[30,10])
from dune.alugrid import aluCubeGrid as leafGridView
grid = leafGridView( cartesianDomain([0,0],[3,1],[30,10]) )

spcU = rtSpace( grid, dimRange=grid.dimension, order=1, storage="fem" )
spcB = bdm( grid, dimRange=grid.dimension, order=1, storage="fem" )

cell  = spcU.cell()
x     = SpatialCoordinate(cell)
exact_u  = as_vector( [x[1] * (1.-x[1]), 0] )

uh = spcU.interpolate( exact_u, name="U" )
ub = spcB.interpolate( exact_u, name="U" )
