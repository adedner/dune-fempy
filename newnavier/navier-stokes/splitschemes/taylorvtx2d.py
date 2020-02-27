
from __future__ import print_function

import math, os, sys, time

import dune.create as create
from dune.grid import cartesianDomain, Marker
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction, replace,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC, expression2GF
import dune.fem
from math import pi,log10
import dune.fem as fem

from dune.fem.function import integrate

from ufl import  sqrt

from dune.fem import parameter

from schemes.scheme import Scheme
from problems.taylorvortex2dclass import TaylorVortex2d


parameter.append({"fem.verboserank": 0, "fem.solver.verbose": 0,
    "fem.solver.errormeassure": "relative",
    "fem.preconditioning" : "true",
    "petsc.kspsolver.method": "gmres",
    "petsc.preconditioning.method": "hypre",
    "istl.preconditioning.method": "ilu",
    "istl.preconditioning.iterations": 0,
    "istl.preconditioning.relaxation": 0.9,
    "istl.gmres.restart": 50,
    "istl.preconditioning.fastilustorage": 1})

method         = Scheme.fromString(parameter.get("method","ipcs"))
useHAdaptivity = parameter.get("hAdaptivity",False)
usePAdaptivity = parameter.get("pAdaptivity",False)
useTAdaptivity = parameter.get("tAdaptivity",False)
maxLevel       = parameter.get("maxLevel",3)
maxOrder       = parameter.get("maxOrder",2)

problem  = TaylorVortex2d

scheme = Scheme(method,problem)

print("method:",Scheme.toString(method))


outputBase = "outputHP" if usePAdaptivity else "outputTaylorvtx"
outputBase += "_"+str(scheme)+str(maxOrder)+str(maxLevel)
# filename   = outputBase+"_"+Scheme.toString(method)
filename = outputBase+"/"+"_"+str(scheme)+str(maxOrder)+str(maxLevel)


# create grid functions from UFL expression of exact solutions
exact_u = expression2GF( problem.grid, problem.exact_u, order=maxOrder+1, name="U")
exact_p = expression2GF( problem.grid, problem.exact_p, order=maxOrder+1, name="P")

solU = problem.solU
solP = problem.solP
spc  = [problem.spcU,problem.spcP]
solution = [solU,solP]
grid = problem.grid
# grid.hierarchicalGrid.globalRefine(4)

timeStep = problem.deltaT
simTime = 0
counter = 0
print("set up output...")
directory = os.path.dirname(filename)
if not os.path.exists(directory):
    os.makedirs(directory)


def polOrder(e,x):
    return problem.spcU.localOrder(e)

vtk = grid.sequencedVTK(filename,
        pointdata={"pressure":solution[1], "exact_p":exact_p},
        celldata={ "order":polOrder},
        pointvector={"velocity":solution[0],"exact_u":exact_u})
vtk()


tolerance = 0.001
gridSize = grid.size(0)
def mark(element):
    solutionLocal = solU.localFunction(element)
    grad = solutionLocal.jacobian(element.geometry.referenceElement.center)
    if grad[0].infinity_norm > 1.2:
        return Marker.refine if element.level < maxLevel else Marker.keep
    else:
        return Marker.coarsen
def markp(element):
    solutionLocal = solU.localFunction(element)
    grad = solutionLocal.jacobian(element.geometry.referenceElement.center)
    if grad[0].infinity_norm > 1.2:
     polorder = 2
    else:
     polorder = 3

    newPolorder = polorder
    return newPolorder


while simTime < problem.endTime:
    print( "Time is:", simTime )
    scheme.solve(simTime,target=solution)
    # grid.hierarchicalGrid.mark(mark)
    # fem.adapt(grid.hierarchicalGrid,[solU,solP])
    # # compute smoothness indicator
    # fem.spaceAdapt(problem.spcU,markp, [solU])
    # fem.loadBalance(grid.hierarchicalGrid,[solU])
    exact_u.t = simTime
    exact_p.t = simTime

    #l2error_fn = dot(solU - exact(simTime)[0], solU - exact(simTime)[0])
    l2error_fn = dot(solU - exact_u, solU - exact_u)
    l2error = sqrt( integrate(grid, l2error_fn, 5) )
    l2error_fn_press = dot(solP - exact_p, solP - exact_p)
    l2errorpress = sqrt( integrate(grid, l2error_fn_press, 5) )
    print('|u_h - u| =', l2error)
    print('|p_h - p| =', l2errorpress)
    simTime += timeStep

    vtk()
