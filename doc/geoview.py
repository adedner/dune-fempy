# <markdowncell>
# <codecell>
import dune.plotting
dune.plotting.block = True
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
from ufl import sqrt, SpatialCoordinate, triangle, as_vector
from dune.grid import structuredGrid
from dune.fem.view import geometryGridView
from dune.fem.function import uflFunction
square = structuredGrid([0,0],[1,1],[10,10])
x = SpatialCoordinate(triangle)
transform = as_vector([ (x[0]+x[1])/sqrt(2), (-x[0]+x[1])*sqrt(2) ])
gridFunction = uflFunction(square,ufl=transform,order=1,name="diamond")
diamond = geometryGridView(gridFunction)
# <markdowncell>
# <codecell>
diamond.plot()
