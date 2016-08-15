from __future__ import print_function

from mpi4py import MPI

import math
import ufl

import dune.models.elliptic
import dune.ufl
import dune.fem

import random

#This is an exact replication of the crystal growth model in dune fem tutorial with one part at a time changed
#used for debuggin

# problem parameters
# ------------------
dimRange     = 2
dimDomain    = 2
maxLevel     = 12
dt           = 5.e-4
endTime      = 0.2
saveinterval = 0.001

#should try and change this so that they don't need ot be made global
global u,v,un


## model taken from www.ctcms.nist.gov/fipy/examples/phase/generated/examples.phase.anisotropy.html
alpha        = 0.015
tau          = 3.e-4
kappa1       = 0.9
kappa2       = 20.
c            = 0.02
N            = 6.


# set up function spaces
uflSpace       = dune.ufl.Space(dimDomain, dimRange)
uflScalarSpace = dune.ufl.Space(dimDomain, 1)
u = ufl.TrialFunction(uflSpace)
v = ufl.TestFunction(uflSpace)
un = ufl.Coefficient(uflSpace)
noise = ufl.Coefficient(uflSpace)

def initial(x):
    r  = (x-[6,6]).two_norm
    return [ 0 if r>0.5 else 1, -0.5 ]
def globalNoise(x):
    return [ 0, 0 ]

#this is the function L(u) in the forcing term
#if unsure return x
def L(x):
    return x + noise[0]

def Sharp(x,y):
    return -2.25*ufl.inner(ufl.grad(x), ufl.grad(y))


# right hand sie (time derivative part + explicit forcing in v)
a_ex = (ufl.inner(un, v) - ufl.inner(un[0], v[1])) * ufl.dx
# left hand side (heat equation in first variable + backward Euler in time)
psi        = ufl.pi/8.0 + ufl.atan_2(ufl.grad(un[0])[1], (ufl.grad(un[0])[0]))
Phi        = ufl.tan(N / 2.0 * psi)
beta       = (1.0 - Phi*Phi) / (1.0 + Phi*Phi)
dbeta_dPhi = -2.0 * N * Phi / (1.0 + Phi*Phi)
fac        = 1.0 + c * beta
diag       = fac * fac
offdiag    = -fac * c * dbeta_dPhi
d0         = ufl.as_vector([diag, offdiag])
d1         = ufl.as_vector([-offdiag, diag])
m          = u[0] - 0.5 - kappa1 / ufl.pi*ufl.atan(kappa2*u[1])
s          = ufl.as_vector([dt / tau * u[0] * (1.0 - u[0]) * m + noise[0], u[0]])

a_im = (alpha*alpha*dt / tau * (ufl.inner(ufl.dot(d0, ufl.grad(u[0])), ufl.grad(v[0])[0]) + ufl.inner(ufl.dot(d1, ufl.grad(u[0])), ufl.grad(v[0])[1]))
       + 2.25 * dt * ufl.inner(ufl.grad(u[1]), ufl.grad(v[1])) + ufl.inner(u,v) - ufl.inner(s,v)) * ufl.dx





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
vtk = grid.writeVTK("crystal", pointdata=[solution], celldata=[grid.levelFunction(), grid.partitionFunction()], number=count)

while t < endTime:
    noise_h.interpolate(noise_gf)
    solution_n.assign(solution)
    scheme.solve(target=solution)
    t += dt
    print('count: ',count,"t = ",t)
    if t > savestep:
        savestep += saveinterval
        count += 1
        vtk.write("crystal", count)
    hgrid.mark(mark)
    hgrid.adapt([solution])
    hgrid.loadBalance([solution])

print("END")
