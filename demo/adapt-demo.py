from __future__ import print_function

# needed on some machines
from mpi4py import MPI

import math

import dune.femmpi as femmpi
import dune.fem as fem

grid = fem.leafGrid("../data/unitcube-2d.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming")

# interpolate some data onto macro grid
phi = grid.interpolate(lambda x: [math.sin(math.pi*x[0])*math.cos(math.pi*x[1])], space="Lagrange", name="phi")

# add phi to vtk output
t = 0
grid.writeVTK( "initial", pointdata=[phi] )

maxLevel = 8
hgrid = grid.hierarchicalGrid
marker = hgrid.marker

def mark(element):
    y = element.geometry.center - [0.5+0.2*math.cos(t), 0.5+0.2*math.sin(t)]
    if y.two_norm2 < 0.2*0.2 and y.two_norm2 > 0.1*0.1:
      return marker.refine if element.level < maxLevel else marker.keep
    else:
      return marker.coarsen

for i in range(0,maxLevel):
    hgrid.mark(mark)
    hgrid.adapt([phi])
    hgrid.loadBalance([phi])

nr = 0
while t < 2*math.pi:
    print('time:',t)
    hgrid.mark(mark)
    hgrid.adapt([phi])
    hgrid.loadBalance([phi])
    phi.interpolate([0])
    grid.writeVTK("adapt", pointdata=[phi], celldata=[grid.levelFunction(), grid.partitionFunction()], number=nr)
    print(femmpi.comm.rank, grid.size(0))
    t = t+0.1
    nr = nr+1
