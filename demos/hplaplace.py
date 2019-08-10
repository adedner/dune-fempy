from __future__ import print_function

import math
from ufl import *

from dune.fem import spaceAdapt, adapt
from dune.grid import cartesianDomain, gridFunction
from dune.alugrid import aluConformGrid as leafGridView
from dune.fem.view import adaptiveLeafGridView as adaptiveGridView
from dune.fem import parameter, doerflerMark, globalRefine
from dune.fem.function import levelFunction
from dune.ufl import Space

import dune.create as create

parameter.append({"fem.verboserank": 0})

domain = cartesianDomain([0,0],[1,1],[16,16])
grid = adaptiveGridView( leafGridView(domain) )
space = create.space("lagrangehp", grid, maxOrder=4, storage="istl")

u = TrialFunction(space)
v = TestFunction(space)
x = SpatialCoordinate(space)
n = FacetNormal(space)
mu = 20 * 16
hT = MaxCellEdgeLength(space)
hS = avg( MaxFacetEdgeLength(space) )
hs = MaxFacetEdgeLength(space)('+')

diffusiveFlux = lambda w,d: d
source = -sin(pi*x[0])*sin(6*pi*x[1])

a  = ( inner(diffusiveFlux(u,grad(u)), grad(v)) + source*v ) * dx
a -= ( inner( outer(jump(u), n('+')), avg(diffusiveFlux(u,grad(v))) ) +\
       inner( avg(diffusiveFlux(u,grad(u))), outer(jump(v), n('+'))) ) * dS
a -= ( inner( outer(u, n), diffusiveFlux(u,grad(v)) ) +\
       inner( diffusiveFlux(u,grad(u)), outer(v, n) ) ) * ds
a += mu/hS * inner(jump(u), jump(v)) * dS
a += mu/hs * inner(u, v) * ds


newtonParameter = {"tolerance": 1e-10, "verbose": "true",
                   "linear.tolerance": 1e-11,
                   "linear.preconditioning.method": "ilu",
                   "linear.preconditioning.iterations": 1, "linear.preconditioning.relaxation": 1.2,
                   "linear.verbose": "false"}
scheme = create.scheme("galerkin", a==0, parameters={"newton." + k: v for k, v in newtonParameter.items()})

solution = space.interpolate([0],name="solution")
scheme.solve(target=solution)
solutionp = solution.copy(name="solutionp")
def markp(element):
    return 1 if element.geometry.center[0]<0.5 else 2
spaceAdapt(space,markp,[solution,solutionp])
grid.writeVTK("pre-hplaplace", pointdata=[solution,solutionp],subsampling=3)
scheme.solve(target=solutionp)
grid.writeVTK("post-hplaplace", pointdata=[solution,solutionp],subsampling=3)

#######################################################################

from dune.fem.space import finiteVolume as estimatorSpace
from dune.fem.operator import galerkin as estimatorOp

fvspace = estimatorSpace(solution.space.grid)
estimate = fvspace.interpolate([0], name="estimate")

u = TrialFunction(space)
v = TestFunction(fvspace)
hT = MaxCellEdgeLength(space)
he = MaxFacetEdgeLength(space)('+')
n = FacetNormal(space)
estimator_ufl = hT**2 *( -div( diffusiveFlux(u,grad(u)) ) + source )**2 * v * dx +\
           he * inner( jump(diffusiveFlux(u,grad(u))), n('+'))**2 * avg(v) * dS +\
         1/he * jump(u)**2 * avg(v) * dS
estimator = estimatorOp(estimator_ufl)
tolerance = 1e-8

@gridFunction(grid,"pDegree")
def pDegree(element,x):
    return space.localOrder(element)

vtk = grid.sequencedVTK("hplaplace", pointdata=[solution,solutionp],
                                     celldata=[estimate,pDegree,levelFunction(grid)],
                                     subsampling=1)
while True:
    estimator(solutionp, estimate)
    eta = math.sqrt( sum(estimate.dofVector) )
    vtk()
    print("estimate:",eta)
    if eta < tolerance:
        break
    marked = doerflerMark(estimate,0.6,layered=0.1)
    # globalRefine(1,solution,solutionp)
    adapt(solution,solutionp)
    spaceAdapt(space,markp,[solution,solutionp])
    print("adapted:",marked)
    scheme.solve(target=solutionp)
