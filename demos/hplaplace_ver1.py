from __future__ import print_function

import math
from ufl import *
import numpy

from dune.fem import spaceAdapt, adapt
from dune.grid import cartesianDomain, gridFunction
from dune.alugrid import aluConformGrid as leafGridView
from dune.fem.view import adaptiveLeafGridView as adaptiveGridView
from dune.fem import parameter, mark, doerflerMark, globalRefine
from dune.fem.function import levelFunction
from dune.ufl import Space, DirichletBC

import dune.create as create

parameter.append({"fem.verboserank": 0})

if False:
    domain = cartesianDomain([0,0],[1,1],[4,4])
else:
    cornerAngle=320
    vertices = numpy.zeros((8, 2))
    vertices[0] = [0, 0]
    for i in range(0, 7):
        vertices[i+1] = [math.cos(cornerAngle/6*math.pi/180*i),
                         math.sin(cornerAngle/6*math.pi/180*i)]
    triangles = numpy.array([[2,1,0], [0,3,2], [4,3,0],
                             [0,5,4], [6,5,0], [0,7,6]])
    domain = {"vertices": vertices, "simplices": triangles}

grid = adaptiveGridView( leafGridView(domain) )
grid.hierarchicalGrid.globalRefine(3)

maxOrder = 4
space   = create.space("lagrangehp", grid, maxOrder=maxOrder, storage="istl")
spacePm = create.space("lagrangehp", grid, maxOrder=maxOrder, storage="istl")

u = TrialFunction(space)
v = TestFunction(space)
x = SpatialCoordinate(space)
n = FacetNormal(space)
mu = 20 * 16
hT = MaxCellEdgeLength(space)
hS = avg( MaxFacetEdgeLength(space) )
hs = MaxFacetEdgeLength(space)('+')

diffusiveFlux = lambda w,d: d
source = 0 # -sin(pi*x[0])*sin(6*pi*x[1])
# exact solution for this angle
Phi = cornerAngle / 180 * pi
phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*pi, 0)
exact = dot(x, x)**(pi/2/Phi) * sin(pi/Phi * phi)

a  = ( inner(diffusiveFlux(u,grad(u)), grad(v)) + source*v ) * dx
a -= ( inner( outer(jump(u), n('+')), avg(diffusiveFlux(u,grad(v))) ) +\
       inner( avg(diffusiveFlux(u,grad(u))), outer(jump(v), n('+'))) ) * dS
a += mu/hS * inner(jump(u), jump(v)) * dS
a -= ( inner( outer(u-exact, n), diffusiveFlux(u,grad(v)) ) +\
       inner( diffusiveFlux(u,grad(u)), outer(v, n) ) ) * ds
a += mu/hs * inner(u-exact, v) * ds

newtonParameter = {"tolerance": 1e-10, "verbose": "true",
                   "linear.tolerance": 1e-11,
                   "linear.preconditioning.method": "ilu",
                   "linear.preconditioning.iterations": 1, "linear.preconditioning.relaxation": 1.2,
                   "linear.verbose": "false"}
scheme = create.scheme("galerkin", [a==0],
          parameters={"newton." + k: v for k, v in newtonParameter.items()})

solution = space.interpolate(0,name="solution")
scheme.solve(target=solution)

#######################################################################

from dune.fem.space import finiteVolume as estimatorSpace
from dune.fem.operator import galerkin as estimatorOp

fvspace = estimatorSpace(grid)
estimate = fvspace.interpolate([0], name="estimate")
estimatePm = fvspace.interpolate([0], name="estimatePm")

u = TrialFunction(space)
v = TestFunction(fvspace)
hT = MaxCellEdgeLength(space)
he = MaxFacetEdgeLength(space)('+')
n = FacetNormal(space)
estimator_ufl = hT**2 *( div( grad(u) ) )**2 * v * dx # +\
# estimator_ufl = hT**2 *( -div( diffusiveFlux(u,grad(u)) ) + source )**2 * v * dx # +\
         #   he * inner( jump(diffusiveFlux(u,grad(u))), n('+'))**2 * avg(v) * dS # +\
         # 1/he * jump(u)**2 * avg(v) * dS +\
         # 1/he * (u-exact)**2 * v * ds
estimator = estimatorOp(estimator_ufl)
tolerance = 1e-8

@gridFunction(grid,"pDegree")
def pDegree(element,x):
    return space.localOrder(element)

solutionPm = spacePm.interpolate(solution)
vtk = grid.sequencedVTK("hplaplace",
     pointdata={"solution":solution,"pm":solutionPm,"exact":exact},
     celldata=[estimate,pDegree,levelFunction(grid)],
     subsampling=1)

####################################################################

pTol = 1e-16
pMarked = [0,]*maxOrder
def markPm(element):
    return max(space.localOrder(element)-1,1)
def markp(element):
    eta = estimate(element,[0,0])[0]
    etaPm = estimatePm(element,[0,0])[0]
    estP = abs(eta-etaPm)
    print(estP)
    polOrder = space.localOrder(element)
    newPolOrder = polOrder
    if estP < pTol:
        newPolOrder = polOrder - 1 if polOrder > 2 else polOrder
    elif estP > 100*pTol:
        newPolOrder = polOrder + 1 if polOrder < maxOrder else polOrder
    pMarked[polOrder-1] += 1
    return polOrder

while True:
    estimator(solution, estimate)
    eta2 = sum(estimate.dofVector)
    vtk()
    print("estimate:",eta2,tolerance)
    if eta2 < tolerance*tolerance*2:
        break
    hMarked = mark(estimate,math.sqrt(eta2)/grid.size(0))
    # globalRefine(1,solution,solution)
    adapt(solution)
    print("h-adapted:",hMarked)

    estimator(solution, estimate)
    spaceAdapt(spacePm,markPm,[solutionPm])
    solutionPm.interpolate(solution)
    estimator(solutionPm, estimatePm)
    spaceAdapt(space,markp,[solution])
    print("p-adapted:",pMarked)

    scheme.solve(target=solution)
