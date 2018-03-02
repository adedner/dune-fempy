from dune.ufl import Space
from ufl import TestFunction, TrialFunction, SpatialCoordinate, ds, dx, inner, grad, dot
import dune.grid
import dune.create as create
from dune.generator import algorithm
from dune.fem.plotting import plotPointData as plot
import dune.common as common
import dune.fem as fem
import math
from dune.fem.function import levelFunction, partitionFunction

def initialRefine(grid,uh,initial_gf):
    hgrid.globalRefine(5)

    for i in range(0,maxLevel):
        hgrid.mark(lambda x: mark(uh,x) )
        fem.adapt(hgrid,[uh])
        uh.interpolate(initial_gf)

def mark(solution,element):

    marker = dune.grid.Marker
    solutionLocal = solution.localFunction(element)
    grad = solutionLocal.jacobian(element.geometry.referenceElement.center)

    if grad.infinity_norm > 0.4  :
        return marker.refine if element.level < maxLevel else marker.keep
    else:
        return marker.coarsen if element.level > minLevel else marker.keep

#First initialise grid using istl
view = create.view("adaptive", grid="ALUConform",constructor = dune.grid.cartesianDomain([-2,-2],[2,2],[3,3]),dimgrid=2)

hgrid = view.hierarchicalGrid

spc = create.space("Lagrange", view, dimrange=1, order=1, storage='istl')

uflSpace = Space(spc)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

#max and mimumal level of grid refinement
maxLevel = 12
minLevel = 5

#setup function with the intial conditions
def initial(x):
    r = (x).two_norm
    ir = 0.5
    return  [ -1 if r > ir else 1]

initial_gf = create.function("global", view, "initial" , 2 , initial)

#uh is the solution at each step of the iteration
uh = spc.interpolate(initial_gf, name="solution")

initialRefine(view,uh,initial_gf)
view.writeVTK("mcf_istl", pointdata=[uh], celldata=[levelFunction(view), partitionFunction(view)], number=0)

#Now initialise grid using fem
spc = create.space("Lagrange", view, dimrange=1, order=1, storage='fem')

uflSpace = Space(spc)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

#uh is the solution at each step of the iteration
uh = spc.interpolate(initial_gf, name="solution")

initialRefine(view,uh,initial_gf)
view.writeVTK("mcf_fem", pointdata=[uh], celldata=[levelFunction(view), partitionFunction(view)], number=0)
