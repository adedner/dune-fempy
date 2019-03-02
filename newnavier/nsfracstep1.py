import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction, replace,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from math import pi,log10
from dune.fem.function import integrate

from ufl import cos, sin, exp, sqrt

from dune.fem import parameter

from burgerclass import Burgers
from stokesclass import Stokes

parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "jacobi", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

Re      = 100.
deltaT  = 0.001
endTime = 0.1
order   = 2

grid = create.grid("ALUCube",constructor=cartesianDomain([-1,-1],[1,1],[50,50]))
spcU = create.space("lagrange", grid, dimrange=grid.dimension, order=order, storage="istl")
spcP = create.space("lagrange", grid, dimrange=1, order=order-1, storage="istl")
x     = SpatialCoordinate(spcU)
time  = NamedConstant(spcU, "t")     # current time

nu    = 1./Re
gamma_t = lambda t: exp(-2.*pi*pi*nu*t)
exact_u = as_vector( [ -1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t(time), sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*nu*time) ] )
exact_p = as_vector( [ -0.25*(cos(2.*pi*x[0])+cos(2.*pi*x[1]))*exp(-4.*pi*pi*nu*time) ] )

exact = lambda t: [replace(exact_u,{time:t}),replace(exact_p,{time:t})]

f           = as_vector( [2.*pi*pi*nu*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t(time) , -2.*pi*pi*nu*sin(1.*pi*x[0])*cos(1.*pi*x[1])*gamma_t(time) ] )
f          += as_vector( [0.25*2.*pi*sin(2.*pi*x[0])*exp(-4.*pi*pi*nu*time),  0.25*2.*pi*sin(2.*pi*x[1])*exp(-4.*pi*pi*nu*time)] )
f          -= as_vector( [nu*2*gamma_t(time)*pi*pi*cos(1.*pi*x[0])*sin(1.*pi*x[1]), nu*-2*gamma_t(time)*pi*pi*sin(1.*pi*x[0])*cos(1.*pi*x[1])] )
f          += as_vector( [cos(1.*pi*x[0])*sin(1.*pi*x[0])*(cos(1.*pi*x[1])*cos(1.*pi*x[1])-sin(1.*pi*x[1])*sin(1.*pi*x[1])),cos(1.*pi*x[1])*sin(1.*pi*x[1])*(cos(1.*pi*x[0])*cos(1.*pi*x[0])-sin(1.*pi*x[0])*sin(1.*pi*x[0]))] )

bcs = [DirichletBC(spcU,[None,None],1)]
##############################################

# TODO: set ufl time constant in expression
solU = spcU.interpolate(exact(0)[0],"velocity")
solP = spcP.interpolate(exact(0)[1],"pressure")
spc  = [spcU,spcP]
solution = [solU,solP]

class PR:
    def __init__(self):
        self.stokes  = Stokes(spc,f,bcs,Re,withBurgers=True)
        self.burgers = Burgers(spc,f,bcs,Re)
    def solve(self,target):
        print( 'Solve step 1 - Stokes' )
        self.stokes.prepare(simTime+0.5*deltaT,2./deltaT,1/2,target)
        self.stokes.solve(target)
        print( 'Solve step 2 - Burgers' )
        self.burgers.prepare(simTime+deltaT,2./deltaT,1/2,target)
        self.burgers.solve(target)
class FractionalTheta:
    def __init__(self):
        self.stokes  = Stokes(spc,f,bcs,Re,withBurgers=True)
        self.burgers = Burgers(spc,f,bcs,Re)
        self.theta = 1-sqrt(2.)/2.
        self.alpha = 0.5
        self.beta  = 1-self.alpha
    def solve(self,target):
        print( 'Solve step 1 - Stokes' )
        self.stokes.prepare(simTime+self.theta*deltaT,1./(self.theta*deltaT),self.alpha,target)
        self.stokes.solve(target)
        print( 'Solve step 2 - Burgers' )
        self.burgers.prepare(simTime+deltaT*(1-self.theta),1/((1.-2*self.theta)*deltaT),self.alpha,target)
        self.burgers.solve(target)
        print( 'Solve step 3 - Stokes' )
        self.stokes.prepare(simTime+deltaT,1./(self.theta*deltaT),self.alpha,target)
        self.stokes.solve(target)

vtk = grid.sequencedVTK("visc1navierfracstep_model2",
        pointdata={"pressure":solution[1], "exact_p":exact(endTime)[1]},
        pointvector={"velocity":solution[0],"exact_u":exact(endTime)[0]})
vtk()

scheme = FractionalTheta() #PR() #
endTime = 0.1
timeStep = deltaT
simTime = 0
counter = 0
while simTime < endTime:
    print( "Time is:", simTime )
    scheme.solve(target=solution)
    simTime += timeStep

    l2error_fn = dot(solU - exact(simTime)[0], solU - exact(simTime)[0])
    l2error = sqrt( integrate(grid, l2error_fn, 5)[0] )
    print('|u_h - u| =', l2error)
    vtk()
