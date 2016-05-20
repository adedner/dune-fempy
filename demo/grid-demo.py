from __future__ import print_function

import math

# dune.fem modules
import dune.fem.grid as grid
import dune.fem.function as function

# just get the grid (only for testing - not used)
onedgrid = grid.leafGrid("../data/unitcube-1d.dgf", "OneDGrid")

for element in onedgrid.elements:
    print( "Center ", element.geometry.center )
    for corner in element.geometry.corners:
        print( "Corner ", corner )

# get the full grid module and then the grid (module needed for grid # functions and output object)
m_yaspgrid = grid.get("YaspGrid", dimgrid=2)
yaspgrid = grid.leafGrid("../data/unitcube-2d.dgf", m_yaspgrid)

vtk_yaspgrid = yaspgrid.vtkWriter()

def expr_global(x):
    return [-(x[1] - 0.5)*math.sin(x[0]*12)]

ggf = yaspgrid.globalGridFunction("expr_global", expr_global)
print("ggf:", ggf, " | ", ggf.name, " with dimRange = ", ggf.dimRange)
ggf.addToVTKWriter(vtk_yaspgrid, vtk_yaspgrid.DataType.PointData)
for element in yaspgrid.elements:
    lf = ggf.localFunction(element)
    x = [0.5, 0.5]
    y = element.geometry.position(x)
    print("ggf( ", y, " ) = ", lf.evaluate(x), " | ", expr_global(y))

def expr_local(element, x):
    geo = element.geometry
    return [abs(expr_global(geo.position(x))[0] - expr_global(geo.center)[0])]

lgf = yaspgrid.localGridFunction("expr_local", expr_local)
lgf.addToVTKWriter(vtk_yaspgrid, vtk_yaspgrid.DataType.CellData)
print("lgf:", lgf, " | ", lgf.name, " with dimRange = ", lgf.dimRange)
for element in yaspgrid.elements:
    lf = lgf.localFunction(element)
    x = [0.5, 0.5]
    y = element.geometry.position(x)
    print("lgf( ", y, " ) = ", lf.evaluate(x), " | ", expr_local(element,x))

ggf = yaspgrid.globalGridFunction("MathExpression", function.MathExpression(["-(x1-0.5)","x0-1./2."]))
ggf.addToVTKWriter(vtk_yaspgrid, vtk_yaspgrid.DataType.PointData)
print("ggf:", ggf, " | ", ggf.name, " with dimRange = ", ggf.dimRange)
for element in yaspgrid.elements:
    lf = ggf.localFunction(element)
    x = [0.5, 0.5]
    y = element.geometry.position(x)
    print("ggf( ", y, " ) = ", lf.evaluate(x))

vtk_yaspgrid.write("grid_demo");
