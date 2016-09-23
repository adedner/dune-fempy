from mpi4py import MPI
from dune.fem import parameter
parameter.append( "../data/parameter" )
parameter.append( {"hallo": 12, "wie": 20, "gehts": "gut?" } )
# print( str(parameter) )
parameter.append( "hallo", 11. )
parameter.append( "mir", "auch" )
parameter.append( "../data/parameter" )
