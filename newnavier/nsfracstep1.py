import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from math import pi,log10
from dune.fem.function import integrate

from ufl import cos, sin, exp, as_vector, dx, grad, inner,sqrt

from dune.fem import parameter

from burgerclass import Burgers
from stokesclass import Stokes

parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

endTime = 0.1
order   = 2

grid = create.grid("ALUCube",constructor=cartesianDomain([-1,-1],[1,1],[50,50]))
spcU = create.space("lagrange", grid, dimrange=grid.dimension, order=order, storage="istl")
spcP = create.space("lagrange", grid, dimrange=1, order=order-1, storage="istl")
x     = SpatialCoordinate(spcU)
dt    = NamedConstant(spcU, "dt")
time  = NamedConstant(spcU, "t")     # current time

gamma_t = lambda t: exp(-2.*pi*pi*mu*t)
exact_u = as_vector( [-1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t(time), sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*mu*time)])
exact_p = as_vector( [ -0.25*(cos(2.*pi*x[0])+cos(2.*pi*x[1]))*exp(-4.*pi*pi*mu*time) ] )

exact_uend = as_vector( [-1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*exp(-2.*pi*pi*mu*0.1), sin(1.*pi*x[0])*cos(1.*pi*x[1])*gamma_t(endTime)

f           = as_vector( [2.*pi*pi*mu*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t(time) , -2.*pi*pi*mu*sin(1.*pi*x[0])*cos(1.*pi*x[1])*gamma_t(time) ] )
f          += as_vector( [0.25*2.*pi*sin(2.*pi*x[0])*exp(-4.*pi*pi*mu*time),  0.25*2.*pi*sin(2.*pi*x[1])*exp(-4.*pi*pi*mu*time)] )
f          -= as_vector( [mu*2*gamma_t(time)*pi*pi*cos(1.*pi*x[0])*sin(1.*pi*x[1]), mu*-2*gamma_t(time)*pi*pi*sin(1.*pi*x[0])*cos(1.*pi*x[1])] )
f          += as_vector( [cos(1.*pi*x[0])*sin(1.*pi*x[0])*(cos(1.*pi*x[1])*cos(1.*pi*x[1])-sin(1.*pi*x[1])*sin(1.*pi*x[1])),cos(1.*pi*x[1])*sin(1.*pi*x[1])*(cos(1.*pi*x[0])*cos(1.*pi*x[0])-sin(1.*pi*x[0])*sin(1.*pi*x[0]))] )
f          += nu*exact_u

##############################################

endTime = 0.1
timeStep = deltaT
time = timeStep
counter = 0
while time < endTime:
    print( "Time is:", time )
    print( 'Solve step 1 - Stokes' )
    solveStokes()
    # print( 'Solve step 2 - Burgers' )
    solveBurger()
    # print( 'Solve step 3 - Stokes' )
    solveStokes()

    l2error_fn = dot(velocity - exact_uend, velocity - exact_uend)
    l2error = sqrt( integrate(grid, l2error_fn, 5)[0] )
    print('|u_h - u| =', l2error)


    vtk()


    time += timeStep

compute()
