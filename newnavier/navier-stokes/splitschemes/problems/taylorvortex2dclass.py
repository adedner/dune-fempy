from ufl import *
import math, os, sys, time

import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction, replace,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from math import pi,log10
from dune.fem.space import lagrange, combined, product

from dune.fem.function import integrate

from ufl import cos, sin, exp, sqrt

from dune.fem import parameter
#
class TaylorVortex2d():
    Re      = 100.
    deltaT  = 0.001
    endTime = 2
    order   =  2
    useClassicalSpaces = True
    useCombSPACE = False
    useDG = False
    useDGPress = False
    useStab = False

    grid = create.view("adaptive", grid="ALUCube",constructor=cartesianDomain([-1,-1],[1,1],[100,100]))
    # spcU = create.space("dgonbhp", grid, dimrange=grid.dimension, order=order, storage="istl")
    if useDG == True:
        spcU = create.space("dgonbhp", grid, dimrange=grid.dimension, order=order, storage="istl")
    else:
        spcU = create.space("lagrange", grid, dimrange=grid.dimension, order=order, storage="istl")

    if useDGPress == True:
        spcP = create.space("dgonbhp", grid, dimrange=1, order=order-1, storage="istl")
    else:
        spcP = create.space("lagrange", grid, dimrange=1, order=order-1, storage="istl")


    spc  = [spcU,spcP]
    x     = SpatialCoordinate(spcU)
    time  = NamedConstant(spcU, "t")     # current time
    nu    = 1./Re
    # gamma_t = lambda t: exp(-2.*pi*pi*1./Re*time)
    exact_u = as_vector( [ -1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*exp(-2.*pi*pi*1./Re*time), sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*nu*time) ] )
    exact_p = as_vector( [ -0.25*(cos(2.*pi*x[0])+cos(2.*pi*x[1]))*exp(-4.*pi*pi*nu*time) ] )
    f       = None
    bcs = [DirichletBC(spcU,[exact_u[0],exact_u[1]],1)]
    bcs_press = [DirichletBC(spcP,[exact_p[0]],1)]
    # solU = spcU.interpolate(exact(0)[0],"velocity")
    # solP = spcP.interpolate(exact(0)[1],"pressure")
    solU = spcU.interpolate([exact_u[0],exact_u[1]],"velocity")
    solP = spcP.interpolate([exact_p[0]],"pressure")
