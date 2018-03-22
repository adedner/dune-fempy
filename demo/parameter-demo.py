from __future__ import print_function
from dune.fem import parameter
parameter.append( "../data/parameter" )
parameter.append( {"hallo": 12, "wie": 20, "gehts": "gut?" } )
# print( str(parameter) )
parameter.append( "hallo", 11. )
parameter["quark"] = 12
parameter.append( "mir", "auch" )
print(parameter["mir"])
try:
    print(parameter["notpresent"])
except KeyError:
    print("notpresent could not be read using __getitem__")
try:
    print(parameter.get("notpresent"))
except KeyError:
    print("notpresent could not be read using get")
print("using get with default value:",int(parameter.get("notpresent",42)))

a = parameter.get("int",10)
bT = parameter.get("boolT",True)
bF = parameter.get("boolF",False)
c = parameter.get("double",2.)
print(type(a),a,type(bT),bT,type(bF),bF,type(c),c)
a = parameter.get("int",-10)
bT = parameter.get("boolT",False)
bF = parameter.get("boolF",True)
c = parameter.get("double",-2.)
print(type(a),a,type(bT),bT,type(bF),bF,type(c),c)

parameter.append( "../data/parameter" )

# test that these changes actually make it through to the C++ side
# problem with determining the C++ type for the string
# from dune.generator import algorithm
# print( algorithm.run('run', 'parameter-demo.hh', "mir", parameter["mir"]) )
# print( algorithm.run('run', 'parameter-demo.hh', "quark", 12 ) )
# print( "not",algorithm.run('run', 'parameter-demo.hh', "quark", -int( parameter["quark"] )) )
