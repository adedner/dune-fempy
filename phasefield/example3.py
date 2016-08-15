from __future__ import print_function

from mpi4py import MPI

import math
import ufl
import ufl.algorithms


import dune.models.elliptic
import dune.ufl
import dune.fem

import random


#these parameters mim the original crystal file
#also has the corret implementation for the interface

# problem parameters
# ------------------
dimRange     = 2
dimDomain    = 2
maxLevel     = 10
dt           = 5.e-4
endTime      = 0.2
saveinterval = 0.001

#should try and change this so that they don't need ot be made global
global u,v,un

# set up function spaces
uflSpace       = dune.ufl.Space(dimDomain, dimRange)
uflScalarSpace = dune.ufl.Space(dimDomain, 1)
u = ufl.TrialFunction(uflSpace)
v = ufl.TestFunction(uflSpace)
un = ufl.Coefficient(uflSpace)
noise = ufl.Coefficient(uflSpace)

def initial(x):
    r  = (x-[6,6]).two_norm
    return [ 0 if r>0.3 else 1, -0.5 ]
def globalNoise(x):
    return [ 0, 0 ]



def setup_phase(eps,K,Tau,sigma,Gamma,L,Sharp):

    #set up the anisotropy symetry
    N = 6

    theta = un[0] # ufl.tan(N / 2.0 * ( ufl.pi/8.0 +  ufl.atan_2(ufl.grad(un[0])[1], ufl.grad(un[0])[0])))
    fac = 1+0.02*ufl.cos(theta) # Gamma(theta)
    tau = Tau(theta)*Gamma(theta)*eps*eps

    theta = ufl.variable(N*theta)
    # facdash = ufl.diff(fac,theta)

    dbeta_dPhi = -N*ufl.sin(N*theta) # -2.0 * N * theta / (1.0 + theta*theta)
    facdash = 0.02*dbeta_dPhi

    #this is the interpolation function
    #can also come in front of K in temperature equations however here using phi_t
    pdash = 6*u[0]*(1-u[0])

    # right hand sie, time derivative part + explicit forcing in v
    a_ex = (ufl.inner(un, v) - K * ufl.inner(un[0], v[1])) * ufl.dx
    # left hand side, heat equation in first variable + backward Euler in time

    #this is so that it can include anisotropic term
    #Made sure as well it is multiplied by Gamma(theta) as the definition of the left hand side is tau(theta)/Gamma(theta)

    #set up the diffusion tensor
    diag       = fac * fac
    offdiag    = -fac * facdash
    d0         = ufl.as_vector([diag, offdiag])
    d1         = ufl.as_vector([-offdiag, diag])



    #including the L(u) term here as well

    #6sqrt2) is the normalisation condtion to make sure \int \phi_{0z}^2 = 6sqrt(2) can be multiplied out at the end of the gibbs thompson condition
    s   = ufl.as_vector([dt / tau * u[0] * (1.0 - u[0]) * (u[0]-0.5) - dt /tau * eps / (6*ufl.sqrt(2)) * L(u[1])* pdash , K * u[0]])


    #I have to change kappa1 = eps*kapp3/6sqrt(2)
    #this is ot normalise it

    #it's easy to add random noise to this term implement this in the future
    a_im = (sigma* eps * eps * dt / tau * (ufl.inner(ufl.dot(d0, ufl.grad(u[0])), ufl.grad(v[0])[0]) + ufl.inner(ufl.dot(d1, ufl.grad(u[0])), ufl.grad(v[0])[1])) - dt * Sharp(u[1],v[1]) + ufl.inner(u,v) - ufl.inner(s,v)) * ufl.dx

    #the end bit here is the sharp interface
    return [a_im, a_ex]



#interface parameter epsilon
#eps = 0.015
eps = 0.015


#This is the latent heat at the interface in our case it is the inverse of the coefficient of u remember this, also in this case the function w(varphi) is varphi and in fact because of the way w is defined this is minus
K = 1.



#this is the coefficient in front of the anisotropic mean curvature term
beta = 1.

alpha = 4./3.

#seems to break when this parameter gets too large
#kappa1 is material parameter in front of u in the gibbs thompson equation
#gamma should be equal to 1000 in this model also need to redefine f though


def Gamma(x):
     #mu here may be way too bit paybe approx 0.02
     return 1.0+0.02* ( (1.0 - x*x) / (1.0 + x*x))

#so the definition of Tau at the moment is tau(theta)*alpha with the alpha taken inside the function
def Tau(x):
    return (1.0/Gamma(x))*alpha

#this is the function L(u) in the forcing term
#if unsure return x
def L(x):
    gamma = 27.00948948
    return gamma*ufl.atan(20.*x)+noise[0]


#define the right hand side of the sharp interface model in variational form and pass this
#here y is the test function
# this corresponds to u_t = Delta u
#can also easily add source term u_t = Delta u + f(u)
def Sharp(x,y):
    return -2.25*ufl.inner(ufl.grad(x), ufl.grad(y))



#returns the implicit and explicit operators
[a_im, a_ex] = setup_phase(eps,K,Tau,beta,Gamma,L,Sharp)


def write():
    from dune.models.elliptic import SourceWriter
    writer = SourceWriter("mymodel.hh")
    writer.openNameSpace('demo')
    model = dune.models.elliptic.compileUFL(a_im == a_ex, tempVars = False)
    model.write(writer, "MyModel")
    writer.closeNameSpace('demo')
    exit(1)
# write()

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
vtk = grid.writeVTK("crystal-151eps", pointdata=[solution], celldata=[grid.levelFunction(), grid.partitionFunction()], number=count)

while t < endTime:
    noise_h.interpolate(noise_gf)
    solution_n.assign(solution)
    scheme.solve(target=solution)
    t += dt
    print('count: ',count,"t = ",t)
    if t > savestep:
        savestep += saveinterval
        count += 1
        vtk.write("crystal-151eps", count)
    hgrid.mark(mark)
    hgrid.adapt([solution])
    hgrid.loadBalance([solution])

print("END")
