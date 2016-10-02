from __future__ import print_function

import math

import dune.fem.function as function

import dune.create as create

# just get the grid (only for testing - not used)
onedgrid = create.grid("OneD", "../data/unitcube-1d.dgf")

for element in onedgrid.elements():
    print( "Center ", element.geometry.center )
    for corner in element.geometry.corners:
        print( "Corner ", corner )

# get the full grid module and then the grid (module needed for grid # functions and output object)
yaspgrid = create.grid("Yasp", "../data/unitcube-2d.dgf", dimgrid=2)
output = []

def expr_global(x):
    return [-(x[1] - 0.5)*math.sin(x[0]*12)]

ggf = create.function("global", yaspgrid, "expr_global", 1, expr_global)
print("ggf:", ggf, " | ", ggf.name, " with dimRange = ", ggf.dimRange)
output.append(ggf)
for element in yaspgrid.elements():
    lf = ggf.localFunction(element)
    x = [0.5, 0.5]
    y = element.geometry.position(x)
    print("ggf( ", y, " ) = ", lf.evaluate(x), " | ", expr_global(y))

def expr_local(element, x):
    geo = element.geometry
    return [abs(expr_global(geo.position(x))[0] - expr_global(geo.center)[0]),
            expr_global(geo.position(x))[0]]

lgf = create.function("local", yaspgrid, "expr_local", 1, expr_local)
output.append(lgf)
print("lgf:", lgf, " | ", lgf.name, " with dimRange = ", lgf.dimRange)
for element in yaspgrid.elements():
    lf = lgf.localFunction(element)
    x = [0.5, 0.5]
    y = element.geometry.position(x)
    print("lgf( ", y, " ) = ", lf.evaluate(x), " | ", expr_local(element,x))

ggf = create.function("global", yaspgrid, "MathExpression", 1,
        lambda x: [-(x[1]-0.5),x[0]-1./2.,x[0],x[1]*x[1],x[0]*x[1],math.sin(x[0]*x[1]),math.exp(-(x[0]-0.5)**2)])
output.append(ggf)
print("ggf:", ggf, " | ", ggf.name, " with dimRange = ", ggf.dimRange)
for element in yaspgrid.elements():
    lf = ggf.localFunction(element)
    x = [0.5, 0.5]
    y = element.geometry.position(x)
    print("ggf( ", y, " ) = ", lf.evaluate(x))

class ExprLocal:
  def __init__(self,gf):
    self.gf_ = gf
  def __call__(self,element, x):
    geo = element.geometry
    return [abs(self.gf_(geo.position(x))[0] - self.gf_(geo.center)[0]),
            self.gf_(geo.position(x))[0]]
lgf = create.function("local", yaspgrid, "ExprLocal", 1, ExprLocal(expr_global))
output.append(lgf)
print("lgf:", lgf, " | ", lgf.name, " with dimRange = ", lgf.dimRange)
for element in yaspgrid.elements():
    lf = lgf.localFunction(element)
    x = [0.5, 0.5]
    y = element.geometry.position(x)
    print("lgf( ", y, " ) = ", lf.evaluate(x), " | ", expr_local(element,x))

# remove all variables in python to make sure they stay alive long enough for the vtk # writer
lgf = "null"
ggf  ="null"

yaspgrid.writeVTK("grid_demo", pointdata=output);
