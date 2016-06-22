from __future__ import print_function

# needed on some machines
from mpi4py import MPI

import math

import dune.fem.grid as grid
import dune.fem.function as function

grid = grid.leafGrid("../data/unitcube-2d.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming")

# interpolate some data onto macro grid
phi = grid.interpolate(lambda x: [math.sin(math.pi*x[0])*math.cos(math.pi*x[1])], space="Lagrange", name="phi", variant="global")

# add phi to vtk output
output = grid.vtkWriter()
phi.addToVTKWriter(output, output.PointData)
output.write("initial")

maxLevel = 8

hgrid = grid.hierarchicalGrid

marker = hgrid.marker
def mark(t,element):
    y = element.geometry.center
    y[0] -= 0.5 + 0.2*math.cos(t)
    y[1] -= 0.5 + 0.2*math.sin(t)
    if y.two_norm2 < 0.1*0.1:
      return marker.refine if element.level < maxLevel else marker.keep
    else:
      return marker.coarsen

t = 0
mark_t = lambda element: mark(t,element)

for i in range(0,maxLevel):
    hgrid.mark(mark_t)
    hgrid.adapt([phi])
    hgrid.loadBalance([phi])

nr = 0
while t < 2*math.pi:
    print('time:',t)
    hgrid.mark(mark_t)
    hgrid.adapt([phi])
    hgrid.loadBalance([phi])
    output.write("adapt"+str(nr))
    print(grid.size(0))
    t = t+0.1
    nr = nr+1
