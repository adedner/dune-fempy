
# coding: utf-8

# # Grid construction and basic usage [(Notebook)][3]
#
# There are a number of ways to construct the macro triangulation of the computational domain. One can either use some avaialble readers using either gmsh files [gmsh][1] or *dune grid format* ([dgf][2]) files. Alternatively, it is possible to directly provide the vertices and grid elements using *numy* arrays.
# Later on we will show some of the most basic method available on the grid ckass.
#
# [1]: http://gmsh.info/doc/texinfo/gmsh.html
# [2]: https://www.dune-project.org/doxygen/2.5.0/group__DuneGridFormatParser.html#details
# [3]: _downloads/grid-construction.ipynb

# In[ ]:


from __future__ import print_function

try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass

try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass

import dune.grid
import dune.create as create


# ### Constructers
# The following shows a simple routing to visualize a dune grid

# In[ ]:


from dune.plotting import block
from matplotlib import pyplot
from matplotlib.collections import PolyCollection
from numpy import amin, amax, linspace

def plotGrid(grid):
    if not grid.dimension == 2:
        print("inline plotting so far only available for 2d grids")
        return
    fig = pyplot.figure()
    polys = grid.polygons()

    for p in polys:   # iterables for different geometry types
        coll = PolyCollection(p,facecolor='none',edgecolor="black",zorder=2)
        pyplot.gca().add_collection(coll)
    fig.gca().set_aspect('equal')
    fig.gca().autoscale()
    pyplot.show(block=block)


# First the simplest approach - this simply results in a cube tesselated with a nicely structured grid:

# In[ ]:


grid = create.grid("ALUCube", dune.grid.cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
plotGrid(grid)
grid = create.grid("ALUConform", dune.grid.cartesianDomain([-1,-1],[1,1],[16,16,16]), dimgrid=2)
plotGrid(grid)


# We now show how describe a simple grid using *numpy* arrays. This approach can then be used to write more complex readers in python to construct grids.

# In[ ]:


import numpy
vertices = numpy.array([(0,0), (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1)])
triangles = numpy.array([(0,1,2), (0,2,3), (0,3,4), (0,4,5), (0,5,6), (0,6,7)])
grid = create.grid("ALUConform", {"vertices": vertices, "simplices": triangles}, dimgrid=2)
plotGrid(grid)


# The next example uses the [Delauny][1] triangulation method availble in *scipy*.
# [1]: https://docs.scipy.org/doc/scipy-0.18.1/reference/tutorial/spatial.html

# In[ ]:


from scipy.spatial import Delaunay

n_angles = 36
n_radii = 8

radii = numpy.linspace(1.0 / n_radii, 1.0, n_radii)
angles = numpy.linspace(0, 2*numpy.pi, n_angles, endpoint=False)
angles = numpy.repeat(angles[..., numpy.newaxis], n_radii, axis=1)

x = numpy.append(0, (radii*numpy.cos(angles)).flatten())
y = numpy.append(0, (radii*numpy.sin(angles)).flatten())

points = numpy.stack((x,y), axis=-1)
triangles = Delaunay(points).simplices

grid = create.grid("ALUConform", {'vertices':points, 'simplices':triangles}, dimgrid=2)
plotGrid(grid)


# Dune provides options for reading grid descriptions from files. First we show how to use the *dune grid format*:
# The grid can be directly described using a python string:

# In[ ]:


dgf = """
INTERVAL
0  0
1  1
16 16
#
"""
grid = create.grid("ALUConform", dune.grid.string2dgf(dgf), dimgrid=2)
plotGrid(grid)


# or providing a [dgf file][1]
# [1]: unitcube-2d.dgf

# In[ ]:


grid = create.grid("ALUCube", (dune.grid.reader.dgf,"unitcube-2d.dgf"), dimgrid=2)
plotGrid(grid)
pyplot.close('all')


# We have just described how to *input* a grid, shown some inline visualization, but we also provide *output* methods for grid (and data) into *vtk* files. For simply writting the grid structure for viewing in e.g. *paraview*:

# In[ ]:


grid.vtkWriter().write("griddemo")


# see [griddemo.vtu][1]
# [1]: griddemo.vtu
#
# ### Basic methods
# Probably the most fundamental information one wants from any grid is the number of entities of different codimension, e,g,, the number of vertices and elements.

# In[ ]:


print("Number of elements:",grid.size(0))
print("Number of vertices:",grid.size(grid.dimension))


# Note that `grid.dimension` give the dimension of the reference elements of the grid while `grid.dimensionworld` returns the dimension of the coordinate space the grid is embedded in.

# In[ ]:


print(grid.dimension,grid.dimensionworld)


# Here is not the place to describe the whole *Dune Grid* interface made available to python. Here is a simple example showing how to iterate over the grid and ascess some simple geometric information for each element:

# In[ ]:


totalArea = 0
for element in grid.elements:
    print( "Center ", element.geometry.center, end=", " )
    print("Corners:", end=' ')
    for corner in element.geometry.corners:
        print(corner, end=' ' )
    print()
    totalArea += element.geometry.volume
print(totalArea)


# Speaking within the Dune context it is important to note that the grid instance `grid` constructor with the methods above is a *leaf grid view*. So even after refinment the grid instance will always provide a view to the finest level of refinement of the hierarchical grid. The full hierarchy grid can be accessed with

# In[ ]:


hgrid = grid.hierarchicalGrid


# The hierarchical grid can be globally refined

# In[ ]:


hgrid.globalRefine(3)
plotGrid(grid)
# is there a global coarsening available?
hgrid.globalRefine(-2)
plotGrid(grid)


# or single elements can be marked for local refinement or coarsening

# In[ ]:


import math
marker = dune.grid.Marker

def mark(element, t):
    y = element.geometry.center - [0.5+0.2*math.cos(t), 0.5+0.2*math.sin(t)]
    if y.two_norm2 < 0.2*0.2 and y.two_norm2 > 0.1*0.1:
      return marker.refine if element.level < 5 else marker.keep
    else:
      return marker.coarsen

for i in range(0,3):
    hgrid.mark(lambda e: mark(e, 0))
    dune.fem.adapt(hgrid)
plotGrid(grid)


# Note that if data is stored on any entity of the grid then it will be valid anymore after calling `globalRefine` or `adapt` on the hierarchical grid as described above. To make it possible to keep for example discrete functions valid during grid modification, these can be passed into the `adapt` method. Note that this requires a special version of the *leaf grid view*:

# In[ ]:


import dune.fem
from dune.fem.space import lagrange
from dune.fem.view import adaptiveLeafGridView
from dune.fem.plotting import plotPointData as plot

adaptiveGV = adaptiveLeafGridView(grid)


# First note that the `plotGrid` function defined above is available in the `dune.fem.plotting` module - with some extra arguments...

# In[ ]:


plot(adaptiveGV)


# In[ ]:


# interpolate some data onto grid
spc = lagrange(adaptiveGV, dimrange=1, order=1)
phi = spc.interpolate(lambda x: [math.sin(math.pi*x[0])*math.cos(math.pi*x[1])], name="phi")

for nr in range(101):
    hgrid.mark(lambda e: mark(e, nr/100.*2.*math.pi))
    dune.fem.adapt(phi)
    if nr % 10 == 0:
        plot(phi)
pyplot.close('all')


# Note the plotting method available in the `dune.fem.plotting` module used in the last example.
#
# ### Special Grid Views
#
# The above required using an `AdaptiveLeafGridView` which allows to efficiently prolong/restrict discrete functions during grid refinement. In this module we provide some other specialized `GridView` which we will now briefly decribe:
# - `FilteredGridView`: a function returning `True` or `False` for each element of grid is passed to the constructor to choose a subset of the elements of the *host grid view* to be part of the new view. The index set associated with this grid view is either (as expected) zero starting and consecutive on the  the new view, or identical to the index set of the original view.
# - `GeometryGridView`: in this case the geometric representation of each element of a *host grid view* can be replaced by providing a [grid function][1] to the constructor.
# [1]: gridfunctions.rst
#
# Lets start with the `FilteredGridView`:

# In[ ]:


from dune.fem.view import filteredGridView
from dune.fem.space import dgonb

# first using a new index set
filter = lambda e: 1 if (e.geometry.center - [0.5, 0.5]).two_norm < 0.25 else 2
subGrid = filteredGridView(grid, filter , 1, True)

# interpolate some data onto the subgrid
spc = lagrange(subGrid, dimrange=1, order=1)
phi = spc.interpolate(lambda x: [math.sin(math.pi*x[0])*math.cos(math.pi*x[1])], name="phi")
plot(phi)

subGrid1 = filteredGridView(grid, filter, 1, True)
subGrid2 = filteredGridView(grid, filter, 2, True)
# interpolate some data onto the subgrid
spc1 = lagrange(subGrid1, dimrange=1, order=1)
spc2 = lagrange(subGrid2, dimrange=1, order=1)
phi1 = spc1.interpolate(lambda x: [-1], name="phi1")
phi2 = spc2.interpolate(lambda x: [ (x-[0.5,0.5] ).two_norm], name="phi2")
spc = dgonb(grid, dimrange=1, order=1)
phi = spc.interpolate(lambda en,x: [phi1.localFunction(en).evaluate(x) if subGrid1.contains(en) else phi2.localFunction(en).evaluate(x)], name="phi")
plot(phi,gridLines="white")


# Note the use of the `contains` method on the view to determin if an element is part of the view or not.
#
# Now a simple example demonstrating hwo to use the `GeometryGridView` to produce  a moving domeain.

# In[ ]:


from dune.fem.view import geometryGridView

t = 0
def expr_global(x):
    return [1.5*x[0],0.5*(x[0]+1.)*x[1]*math.cos(0.1+2.*math.pi*t)]

gf = create.function("global", grid, "coordinates", 1, expr_global)
spc = create.space("lagrange", grid, dimrange=2, order=1)
df = spc.interpolate(gf, name="test")

geogrid = geometryGridView(df)
gfnew = create.function("global", geogrid, "expression", 1, expr_global)

dt = 0.01
count = 0
while t < 1:
    t += dt
    count += 1
    df.interpolate(gf)
    if count%10 == 0:
        plot(gfnew,xlim=[-0.1,1.6],ylim=[-1.25,1.1])
# now one plot without the grid
plot(gfnew,gridLines="")
pyplot.close('all')
