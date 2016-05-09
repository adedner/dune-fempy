from __future__ import print_function
import math

import dune.fem.grid as grid
import dune.fem.gridfunction as gridfunction

grid = grid.leafGrid("../data/unitcube-2d.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming")

# interpolate some data onto macro grid
#phi = grid.interpolate(grid.getGlobal("phi", gridfunction.MathExpression(["math.sin(math.pi*x0)*math.cos(math.pi*x1)"])), "phi", polorder=2)

# add phi to vtk output
output = grid.vtkOutput()
#output.add(phi)
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

hgrid.loadBalance()
for i in range(0,maxLevel):
    hgrid.mark(mark_t)
    hgrid.adapt()
    hgrid.loadBalance()
    #hgrid.adapt([phi])
    #hgrid.loadBalance([phi])

nr = 0
while t < 2*math.pi:
    print('time:',t)
    hgrid.mark(mark_t)
    hgrid.adapt()
    hgrid.loadBalance()
    #hgrid.adapt([phi])
    #hgrid.loadBalance([phi])
    output.write("adapt"+str(nr))
    print(grid.size(0))
    t = t+0.1
    nr = nr+1
