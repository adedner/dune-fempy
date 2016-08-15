from __future__ import print_function

from mpi4py import MPI

import math
import ufl
from ufl import *

import dune.models.elliptic
import dune.ufl
import dune.fem

import random

#this is needed for definiing exp in the numerical integration
import numpy as np

#for integration of function v

from scipy.integrate import quad

#This file is for the isotropic case and the convergence testing
#model taken from 6.4
#http://www.sciencedirect.com/science/article/pii/S0021999196900959
#Computation of Three Dimensional Dendrites with Finite Elements


# problem parameters
# ------------------
dimRange     = 2
dimDomain    = 2
maxLevel     = 12
dt           = 5.e-4
endTime      = 1.0
saveinterval = 0.001

#time constant to keep track of time
time_const = Constant(triangle)

#should try and change this so that they don't need ot be made global
global u,v,un

# set up function spaces
uflSpace       = dune.ufl.Space(dimDomain, dimRange)
uflScalarSpace = dune.ufl.Space(dimDomain, 1)
u = ufl.TrialFunction(uflSpace)
v = ufl.TestFunction(uflSpace)
un = ufl.Coefficient(uflSpace)
noise = ufl.Coefficient(uflSpace)

#setting the time step as a constant locally on a triangle so that it can be changed at each timestep

def initial(x):
    z  = (x-[6,6]).two_norm
    #return [ 0 if z>0.5 else 1, -0.5]
    #new initial condition of the exact solution]
    return [ 0 if z> r_exact(0)  else 1, w_exact(0) if z < r_exact(0) else w_exact(0) + v_exact(z/r_exact(0))  ]

def globalNoise(x):
    return [ 0, 0 ]

#takes a numer usually the radius divded by r(t) and returns expression
def v_exact(s):
    #define the integral
    def integrand(x):
        return  np.exp(-0.25*x*x)/x

    I = -np.exp(0.25)/2.*quad(integrand ,1,s)[0]
    #return the first arugment of the integral which is the result
    #the second argument is the estimated error bound
    return I


#defines the function w for the exact solution also used in initial conditions
def w_exact(t):
    return -0.02*3/(2*r_exact(t))

#gives the exact radius of the circle at each time
def r_exact(t):
    return math.sqrt(0.5*0.5+t)


def setup_phase(eps,K,Tau,sigma,Gamma,L,Sharp):

    theta      = ufl.atan_2(ufl.grad(un[0])[1], (ufl.grad(un[0])[0]))

    #this is the interpolation function
    #can also come in front of K in temperature equations however here using phi_t
    pdash = 6*u[0]*(1-u[0])

    # right hand sie, time derivative part + explicit forcing in v
    a_ex = (ufl.inner(un, v) - K * ufl.inner(un[0], v[1])) * ufl.dx
    # left hand side, heat equation in first variable + backward Euler in time

    fac = Gamma(theta)
    fac = 1.



    #this is so that it can include anisotropic term
    #Made sure as well it is multiplied by Gamma(theta) as the definition of the left hand side is tau(theta)/Gamma(theta)
    tau = Tau(theta)*Gamma(theta)*eps*eps

    #theta = ufl.variable(theta)

    #facdash = ufl.diff(fac,theta)
    facdash = 0.

    #set up the diffusion tensor
    diag       = fac * fac
    offdiag    = -fac * facdash
    d0         = ufl.as_vector([diag, offdiag])
    d1         = ufl.as_vector([-offdiag, diag])



    #including the L(u) term here as well

    #6sqrt2) is the normalisation condtion to make sure \int \phi_{0z}^2 = 6sqrt(2) can be multiplied out at the end of the gibbs thompson condition
    s   = ufl.as_vector([dt / tau * u[0] * (1.0 - u[0]) * (u[0]-0.5) - dt /tau *  eps / (6*ufl.sqrt(2)) * L(u[1])* pdash , K * u[0]])


    #I have to change kappa1 = eps*kapp3/6sqrt(2)
    #this is ot normalise it

    #it's easy to add random noise to this term implement this in the future
    a_im = (sigma* eps * eps * dt / tau * (ufl.inner(ufl.dot(d0, ufl.grad(u[0])), ufl.grad(v[0])[0]) + ufl.inner(ufl.dot(d1, ufl.grad(u[0])), ufl.grad(v[0])[1])) - dt * Sharp(u[1],v[1]) + ufl.inner(u,v) - ufl.inner(s,v)) * ufl.dx

    #the end bit here is the sharp interface

    return [a_im, a_ex]


#interface parameter epsilon
eps = 0.005

#This is the latent heat at the interface in our case it is the inverse of the coefficient of u remember this, also in this case the function w(varphi) is varphi and in fact because of the way w is defined this is minus
K = 1.


#this is the coefficient in front of the anisotropic mean curvature term
beta = 1.

#seems to break when this parameter gets too large
#kappa1 is material parameter in front of u in the gibbs thompson equation
#gamma should be equal to 1000 in this model also need to redefine f though


#this is the anisotropic mean curvature function
def Gamma(x):
    return 1.0

#so the definition of Tau at the moment is tau(theta)*alpha with the alpha taken inside the function
def Tau(x):
    alpha = 1.0
    #return 1.0/Gamma(x)
    return 1.0 * alpha

#this is the function L(u) in the forcing term
#if unsure return x
def L(x):
    gamma = 50.
    return gamma*x + noise[0]

#define the right hand side of the sharp interface model in variational form and pass this
#here y is the test function
# this corresponds to u_t = Delta u
#can also easily add source term u_t = Delta u + f(u)
def Sharp(x,y):
    f = ((0.02*3)/4)*ufl.algebra.Power(0.5*0.5+time_const,-3./2.)

    #should there be a dx term at the end of this surely it wouldn't work otherwise?
    return -ufl.inner(ufl.grad(x), ufl.grad(y)) + ufl.inner(f,y)

#returns the implicit and explicit operators
[a_im, a_ex] = setup_phase(eps,K,Tau,beta,Gamma,L,Sharp)


# basic setup
# -----------
grid       = dune.fem.leafGrid("../data/crystal-2d.dgf", "ALUSimplexGrid", dimgrid=dimDomain, refinement="conforming")
spc        = dune.fem.create.space("Lagrange", grid, dimrange=dimRange, polorder=1)
initial_gf = grid.globalGridFunction("initial", initial)
noise_gf   = grid.globalGridFunction("noise", globalNoise)
solution   = spc.interpolate(initial_gf, name="solution")
solution_n = spc.interpolate(initial_gf, name="solution_n")
noise_h    = spc.interpolate(noise_gf, name="noise_n")

# setup scheme
# ------------
model  = dune.models.elliptic.importModel(grid, dune.models.elliptic.compileUFL(a_im == a_ex)).get()
scheme = dune.fem.create.scheme("FemScheme", solution, model, "scheme")

print(noise.count())
model.setCoefficient(un, solution_n)
model.setCoefficient(noise, noise_h)

# marking strategy
# ----------------
def mark(element):
    solutionLocal = solution.localFunction(element)
    grad = solutionLocal.jacobian(element.geometry.domain.center)
    if grad[0].infinity_norm > 1.0:
      return hgrid.marker.refine if element.level < maxLevel else hgrid.marker.keep
    else:
      return hgrid.marker.coarsen

# initial grid refinement
# -----------------------
hgrid = grid.hierarchicalGrid
grid.globalRefine(2)
for i in range(0,maxLevel):
    hgrid.mark(mark)
    hgrid.adapt([solution])
    hgrid.loadBalance([solution])
    solution.interpolate(initial_gf)

# time loop
# ---------
count    = 0
t        = 0.0
savestep = saveinterval

#simulation name for all the files
sim_name = 'crystal_exact_05eps'

vtk = grid.writeVTK(sim_name, pointdata=[solution], celldata=[grid.levelFunction(), grid.partitionFunction()], number=count)

#open file for writing
file = open(sim_name, 'w')


#set the initial value of time
model.setConstant(time_const, [0.])

#the problem is time_const is a dune constant type we need it as a float

def exact(x):
    #t = time_const[0]
    z  = (x-[6,6]).two_norm
    #the first argument being the phase field doesn't matter
    return [ 0 , w_exact(t) if z < r_exact(t) else w_exact(t) + v_exact(z/r_exact(t))  ]


# function used for computing approximation error
def l2error(en,x):
    #what was solution here was initially uh
    y = en.geometry.position(x)
    val = solution.localFunction(en).evaluate(x) - exact(y)
    return [ math.sqrt( val[1]*val[1]) ];

l2error_gf = grid.localGridFunction( "error", l2error )


while t < endTime:
    noise_h.interpolate(noise_gf)
    solution_n.assign(solution)
    scheme.solve(target=solution)
    t += dt

    #set the new value of time in the code so that this parameter can be included in the source term
    #not quite sure why it's inside of the brackets but seems to be working?
    model.setConstant(time_const, [t])


    error = grid.l2Norm(l2error_gf)
    print('count: ',count,"t = ",t,'error=',error)
    #outputstring = ('t = ',t,'error=',error, '\n')
    file.write(str(t) + ',' + str(error) + '\n')

    if t > savestep:
        savestep += saveinterval
        count += 1
        vtk.write(sim_name, count)
    hgrid.mark(mark)
    hgrid.adapt([solution])
    hgrid.loadBalance([solution])

print("END")
