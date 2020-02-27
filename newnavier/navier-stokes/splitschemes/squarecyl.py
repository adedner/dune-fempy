
from __future__ import print_function

import math, os, sys, time

import dune.create as create
from dune.grid import cartesianDomain, Marker
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction, replace,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC, expression2GF
import dune.fem as fem
from math import pi,log10

from dune.fem.function import integrate

from ufl import cos, sin, exp, sqrt

from dune.fem import parameter

from schemes.scheme import Scheme
from problems.squarecylinderclass import SquareCylinder


parameter.append({"fem.verboserank": 0, "fem.solver.verbose": 0,
    "fem.solver.errormeassure": "relative",
    "fem.preconditioning" : "true",
    "petsc.kspsolver.method": "gmres",
    "petsc.preconditioning.method": "hypre",
    "istl.preconditioning.method": "amg-ilu",
    "istl.preconditioning.iterations": 0,
    "istl.preconditioning.relaxation": 0.9,
    "istl.gmres.restart": 50,
    "istl.preconditioning.fastilustorage": 1})

method         = Scheme.fromString(parameter.get("method","ipcs"))
useHAdaptivity = parameter.get("hAdaptivity",False)
usePAdaptivity = parameter.get("pAdaptivity",False)
useTAdaptivity = parameter.get("tAdaptivity",False)
maxLevel       = parameter.get("maxLevel",5)
maxOrder       = parameter.get("maxOrder",3)

problem  = SquareCylinder

scheme = Scheme(method,problem)


print("method:",Scheme.toString(method))



outputBase = "outputHP" if usePAdaptivity else "outputKarmanvtx"
outputBase += "_"+str(scheme)+str(maxOrder)+str(maxLevel)
# filename   = outputBase+"_"+Scheme.toString(method)
filename = outputBase+"/"+"_"+str(scheme)+str(maxOrder)+str(maxLevel)


solU = problem.solU
solP = problem.solP
spc  = [problem.spcU,problem.spcP]
solution = [solU,solP]
grid = problem.grid
grid.hierarchicalGrid.globalRefine(2)

timeStep = problem.deltaT
simTime = 0
counter = 0
print("set up output...")
# directory = os.path.dirname(filename)
# if not os.path.exists(directory):
#     os.makedirs(directory)

# vtk = grid.sequencedVTK(filename,
#         pointdata={"pressure":solution[1]},
#         pointvector={"velocity":solution[0]})
# vtk()


# def polOrder(e,x):
#     return problem.spcU.localOrder(e)

vtk = grid.sequencedVTK('parall',
        pointdata={"pressure":solution[1]},
        # celldata={ "order":polOrder},
        pointvector={"velocity":solution[0]})
# vtk()


# tolerance = 0.001
# gridSize = grid.size(0)
# def mark(element):
#     solutionLocal = solU.localFunction(element)
#     grad = solutionLocal.jacobian(element.geometry.referenceElement.center)
#     if grad[0].infinity_norm > 1.2:
#         return Marker.refine if element.level < maxLevel else Marker.keep
#     else:
#         return Marker.coarsen
# def markp(element):
#     # estimateLocal   = estimate.localFunction(element)
#     # estimateLocalp1 = estimate_pm1.localFunction(element)
#     solutionLocal = solU.localFunction(element)
#     grad = solutionLocal.jacobian(element.geometry.referenceElement.center)
#     if grad[0].infinity_norm > 1.2:
#      polorder = 2
#     else:
#      polorder = 3
#     # get current polorder
#     # evaluate smoothness indicator
#     # eta = estimateLocal.evaluate(element.geometry.domain.center) / estimateLocalp1.evaluate(element.geometry.domain.center)
#     # newPolorder = polorder
#     # if eta[0] > ptol:
#     #   newPolorder = polorder-1 if polorder > 1 else polorder
#     # else:
#     #   newPolOrder = polorder+1 if polorder < spc.order else polorder
#     # # store new polorder in estimate_pm1
#     # estimateLocalp1[ 0 ] = newPolorder
#     newPolorder = polorder
#     return newPolorder


while simTime < problem.endTime:
    print( "Time is:", simTime )
    scheme.solve(simTime,target=solution)
    # grid.hierarchicalGrid.mark(mark)
    # fem.adapt(grid.hierarchicalGrid,[solU,solP])
    # # compute smoothness indicator
    # fem.spaceAdapt(problem.spcU,markp, [solU])
    # fem.loadBalance(grid.hierarchicalGrid,[solU])

    # l2error_fn = dot(solU - exact(simTime)[0], solU - exact(simTime)[0])
    # l2error = sqrt( integrate(grid, l2error_fn, 5)[0] )
    # print('|u_h - u| =', l2error)
    vtk()
    simTime += timeStep

    # grid.writeVTK("dg-laplace", celldata={"velocity":solution[0]})
