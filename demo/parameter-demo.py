from __future__ import print_function
from dune.fem import parameter
parameter.append( "../data/parameter" )
parameter.append( {"hallo": 12, "wie": 20, "gehts": "gut?" } )
# print( str(parameter) )
parameter.append( "hallo", 11. )
parameter["quark"] = 12
parameter.append( "mir", "auch" )
print(parameter["mir"])
parameter.append( "../data/parameter" )

# test that these changes actually make it through to the C++ side
from dune.generator import algorithm
print( algorithm.run('run', 'parameter-demo.hh', unicode("mir"), parameter["mir"]) )
print( algorithm.run('run', 'parameter-demo.hh', "quark", 12 ) )
print( "not",algorithm.run('run', 'parameter-demo.hh', "quark", -int( parameter["quark"] )) )
