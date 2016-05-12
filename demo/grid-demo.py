from __future__ import print_function

import math

# dune.fem modules
import dune.fem.grid as grid
#import dune.fem.gridfunction as gf

# just get the grid (only for testing - not used)
onedgrid = grid.leafGrid("../data/unitcube-1d.dgf", "OneDGrid")

for element in onedgrid.elements:
    print( "Center ", element.geometry.center )
    for corner in element.geometry.corners:
        print( "Corner ", corner )

# get the full grid module and then the grid (module needed for grid # functions and output object)
m_yaspgrid = grid.get("YaspGrid", dimgrid=2)
yaspgrid = m_yaspgrid.LeafGrid("../data/unitcube-2d.dgf")

def expr_global(x):
    return [-(x[1] - 0.5)*math.sin(x[0]*12)]

ggf = yaspgrid.globalGridFunction("expr_global", expr_global)
print("ggf:", ggf.name, " with dimRange = ", ggf.dimRange)
for element in yaspgrid.elements:
    lf = ggf.localFunction(element)
    x = [0.5, 0.5]
    y = element.geometry.position(x)
    print("ggf( ", y, " ) = ", lf.evaluate(x), " | ", expr_global(y))

def expr_local(element, x):
  geo = element.geometry
  return [abs(expr_global(geo.position(x))[0] - expr_global(geo.center)[0])]

lgf = yaspgrid.localGridFunction("expr_local", expr_local)
print("lgf:", lgf.name, " with dimRange = ", lgf.dimRange)
for element in yaspgrid.elements:
    lf = lgf.localFunction(element)
    x = [0.5, 0.5]
    y = element.geometry.position(x)
    print("lgf( ", y, " ) = ", lf.evaluate(x), " | ", expr_local(element,x))

# ... a 3d alugrid
#start_time = timeit.default_timer()
#grid3 = grid.leafGrid("../data/unitcube-3d.dgf","ALUSimplexGrid", dimgrid=3, refinement="conforming")
#print("Building ALUGrid took ", timeit.default_timer() - start_time, "s")

# refine a bit
#grid2.globalRefine(4)
#grid3.globalRefine(6)

# define two global functions
#def expr_global(x,r):
#    r[0] = -(x[1]-0.5)*math.sin(x[0]*12)
#expr1 = gf.FuncExpression(1,expr_global)
#expr2 = gf.MathExpression(["-(x1-0.5)","x0-1./2."])

# define one local function
#def expr_local(en,x,r):
#  y = en.geometry().position(x)
#  expr_global(y,r)
#  m = en.geometry().position([1./3.,1./3.])
#  rm = [0]
#  expr_global(m,rm)
#  r[0] = abs(r[0]-rm[0])
#expr3 = gf.LocalFuncExpression(1,expr_local)

# convert the global expressions to grid functions on grid2
#scalarGlobal = alugrid2d.getGlobal("scalar",expr1)
#vectorGlobal = alugrid2d.getGlobal("vector",expr2)
#scalarLocal  = alugrid2d.getLocal("local",expr3)
# get an output object from the 2d grid module and write vtk file
#out2 = alugrid2d.VTKOutput( grid2 )
#out2.add( scalarGlobal )
#out2.add( vectorGlobal )
#out2.add( scalarLocal )
#out2.write( "out2d" )

# convert the global expressions to grid functions on grid3
#scalarGlobal = grid3.getGlobal("scalar",expr1)
#vectorGlobal = grid3.getGlobal("vector",expr2)
#scalarLocal  = grid3.getLocal("local",expr3)
# get an output object from the 3d grid module and write vtk file
#out3 = grid3.vtkOutput()
#out3.add( scalarGlobal )
#out3.add( vectorGlobal )
#out3.add( scalarLocal )  # this is actually a mistake since m has only 2 entries
#out3.write( "out3d" )

#out3a = grid3.vtkOutput()
#out3a.write( "grid3d" )

# need to catch
# out2.add( scalarGlobal )
