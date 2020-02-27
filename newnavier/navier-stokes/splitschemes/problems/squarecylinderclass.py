from ufl import *
import math, os, sys, time
from dune.grid import reader

import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction, replace,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from math import pi,log10
from dune.fem.space import lagrange, combined, product

from dune.fem.function import integrate
from dune.alugrid import aluCubeGrid as leafGridView

from ufl import cos, sin, exp, sqrt

from dune.fem import parameter
#

class SquareCylinder():
    Re      = 1000.
    nu    = 1./Re
    deltaT  = 0.00075
    # deltaT  = 0.001
    endTime = 5
    order   = 2
    useClassicalSpaces = True
    useCombSPACE = False
    useDG = False
    useDGPress = False
    useStab = False
    # grid = create.view("adaptive", grid="ALUCube", constructor="../../data/channelflow.dgf", dimgrid=2)

    # grid = create.grid("ALUCube",constructor=cartesianDomain([-1,-1],[1,1],[50,50]))
    # grid = create.view("adaptive", grid="ALUCube", constructor="../../data/karmanvortexstreet_1.msh", dimgrid=2)
    domain = (reader.gmsh, "../../data/cylinder.msh")
    grid  = leafGridView( domain, dimgrid=2 )


    if useCombSPACE == True and useClassicalSpaces == False:
        if useDG == True:
            spc  = combined( dgspcU, dgspcP )
        else:
            spc  = combined( spcU, spcP )
    elif useCombSPACE == False and useClassicalSpaces == False:
        if useDG == True:
            spc  = create.space("dgonbhp", grid, dimrange=3, order=1, storage="istl")
        else:
            spc  = create.space("lagrange", grid, dimrange=3, order=2, storage="istl")
    else:
        if useDG == True:
            spcU = create.space("dgonbhp", grid, dimrange=2, order=order, storage="fem")
            spcP = create.space("lagrange", grid, dimrange=1, order=order-1, storage="fem")

            spc  = [spcU,spcP]
            x     = SpatialCoordinate(spcU)
            time  = NamedConstant(spcU, "t")     # current time

            # gamma_t = lambda t: exp(-2.*pi*pi*1./Re*t)
            exact_u = as_vector([1.5/(0.41*0.41)*(0.41-x[1])*4*x[1]* conditional(x[0]<1e-3,1.,0.),0])
            exact_p = as_vector( [ 0 ] )
            f           = None
            solU = spcU.interpolate([exact_u[0],exact_u[1]],"velocity")
            solP = spcP.interpolate([exact_p[0]],"pressure")
            bcs = [DirichletBC(spcU, [None,None], 2),       # bottom/top
                   DirichletBC(spcU, [None,None], 4),       # bottom/top
                   DirichletBC(spcU, [1.5/(0.41*0.41)*(0.41-x[1])*4*x[1],0], 3),             # left
                   DirichletBC(spcU, [0,0], 1)]
            bcs_press = [DirichletBC(spcP, [0], 2),       # bottom/top
                         DirichletBC(spcP, [0], 4),       # bottom/top
                         DirichletBC(spcP, [0], 3),             # left
                         DirichletBC(spcP, [0], 1)]
        else:
            spcU = create.space("lagrange", grid, dimrange=2, order=order, storage="fem")
            spcP = create.space("lagrange", grid, dimrange=1, order=order-1, storage="fem")
            spc  = [spcU,spcP]

            x     = SpatialCoordinate(spcU)
            time  = NamedConstant(spcU, "t")     # current time

            # 4*1.5*x[1]*(0.41-x[1]*sin(pi*time*1/8)*1/(0.41*0.41*0.41*0.41))
            # gamma_t = lambda t: exp(-2.*pi*pi*1./Re*t)
            exact_u = as_vector([1.5/(0.41*0.41)*(0.41-x[1])*4*x[1]* conditional(x[0]<1e-3,1.,0.),0])
            # exact_u = as_vector([4*1.5*x[1]*(0.41-x[1])*sin(pi*time*1/8)*1/(0.41*0.41*0.41*0.41)* conditional(x[0]<1e-3,1.,0.),0])
            exact_p = as_vector( [ 0 ] )
            f           = None
            bcs = [DirichletBC(spcU, exact_u, conditional(x[0]<2.19,1.,0.)),
                   DirichletBC(spcU, [None,None],conditional(x[0]>2.19,1.,0.))
                #    DirichletBC(spcU, [0,0], conditional(x[1]>0.40,1.,0.)),
                #    DirichletBC(spcU, [0,0], conditional(x[1]<1e-3,1.,0.)),
                #    DirichletBC(spcU, [1.5/(0.41*0.41)*(0.41-x[1])*4*x[1],0], conditional(x[0]<1e-3,1.,0.)),             # left
                   #    DirichletBC(spcU, [4*1.5*x[1]*(0.41-x[1])*sin(pi*time*1/8)*1/(0.41*0.41*0.41*0.41),0], 3),             # left
                   ]

            bcs_press = [DirichletBC(spcP, [None], 2),       # bottom/top
                         DirichletBC(spcP, [None], 5),
                         DirichletBC(spcP, [0], conditional(x[0]>2.19,1.,0.)),       # bottom/top
                         DirichletBC(spcP, [None], 3),             # left
                         DirichletBC(spcP, [None], 6),
                         DirichletBC(spcP, [None], 7),       # bottom/top
                         DirichletBC(spcP, [None], 8),             # left
                         DirichletBC(spcP, [None], 9)]

            outflowVeldirich = conditional(x[0]>7.999,0.,1. )

            dotdirich = conditional(x[0]+1<1e-3,1.,0. )
            pressdirich = conditional(x[0]>4.999,1.,0. )

            solU = spcU.interpolate([exact_u[0],exact_u[1]],"velocity")
            solP = spcP.interpolate([exact_p[0]],"pressure")
