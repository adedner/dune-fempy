from __future__ import print_function

import math

import dune.fem

from dune.grid import cartesianDomain, Marker
from dune.alugrid import aluConformGrid
from dune.fem.function import levelFunction, partitionFunction
from dune.fem.space import lagrange
from dune.fem.view import adaptiveLeafGridView

domain = cartesianDomain([0, 0], [1, 1], [8, 8])
grid = adaptiveLeafGridView(aluConformGrid(domain, dimgrid=2))

# interpolate some data onto macro grid
spc = lagrange(grid, dimrange=1, order=1)
phi = spc.interpolate(lambda x: [math.sin(math.pi*x[0])*math.cos(math.pi*x[1])], name="phi")

# add phi to vtk output
grid.writeVTK("initial", pointdata=[phi])

maxLevel = 8
hgrid = grid.hierarchicalGrid
marker = Marker

def mark(element, t):
    y = element.geometry.center - [0.5+0.2*math.cos(t), 0.5+0.2*math.sin(t)]
    if y.two_norm2 < 0.2*0.2 and y.two_norm2 > 0.1*0.1:
      return marker.refine if element.level < maxLevel else marker.keep
    else:
      return marker.coarsen

for i in range(0,maxLevel):
    hgrid.mark(lambda e: mark(e, 0))
    dune.fem.adapt(hgrid, [phi])
    dune.fem.loadBalance(hgrid, [phi])

vtk = grid.sequencedVTK("adapt", pointdata=[phi], celldata=[levelFunction(grid,name="gridLevel"), partitionFunction(grid)])
vtk()
t = 0
while t < 2*math.pi:
    if grid.comm.rank == 0:
        print('time:', t)
    hgrid.mark(lambda e: mark(e, t))
    dune.fem.adapt(hgrid, [phi])
    dune.fem.loadBalance(hgrid, [phi])
    vtk()
    print("[" + str(grid.comm.rank), "] Size: " + str(grid.size(0)))
    t += 0.1
