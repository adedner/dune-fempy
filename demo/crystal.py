from __future__ import print_function
import math
from mpi4py import MPI

import ufl
import dune.models.femufl as duneuflmodel
import dune.fem as fem
import dune.fem.gridpart as gridpart
import dune.fem.space as space
import dune.fem.scheme as scheme
import dune.fem.function as function

#################################################################
dimRange     = 2
dimDomain    = 2
maxLevel     = 10
dt           = 5.e-4
endTime      = 0.1
saveinterval = 0.001

## model taken from www.ctcms.nist.gov/fipy/examples/phase/generated/examples.phase.anisotropy.html
alpha        = 0.015
tau          = 3.e-4
kappa1       = 0.9
kappa2       = 20.
c            = 0.02
N            = 6.
def initial(x):
    r  = (x-[6,6]).two_norm
    return [ 0 if r>0.1 else 1, -0.5 ]

#################################################################
# Basic setup
# -----------
# set up reference domain
grid2d    = fem.leafGrid("../data/crystal-2d.dgf", "ALUSimplexGrid", dimgrid=dimDomain, refinement="conforming")
grid2d.hierarchicalGrid.globalRefine(3)
sp        = space.create("Lagrange", grid2d, dimrange=dimRange, polorder=1)
level_gf  = grid2d.localGridFunction("level", function.Levels())

# set up left and right hand side models
# --------------------------------------
model     = duneuflmodel.DuneUFLModel(2,dimDomain)
u         = model.trialFunction()
v         = model.testFunction()

# right hand sie (time derivative part + explicit forcing in v)
a = ( ufl.inner(u,v) - ufl.inner(u[0],v[1]) ) * ufl.dx(0)
model.generate(a)
# now add implicit part of forcing to source
rhsModel = model.makeAndImport(grid2d,name="crystal_right").get()

model.clear()

# left hand side (heat equation in first variable + backward Euler in time)
un         = model.coefficient('un',dimRange)
psi        = ufl.pi/8. + ufl.atan( ufl.grad(un[0])[1] / (ufl.grad(un[0])[0]+1e-8) )
Phi        = ufl.tan(N/2.*psi)
beta       = ( 1. - Phi*Phi ) / (1. + Phi*Phi)
dbeta_dPhi = -2. * N * Phi / (1. + Phi*Phi)
fac        = 1. + c*beta
diag       = fac*fac
offdiag    = -fac * c * dbeta_dPhi
d0         = ufl.as_vector([diag, offdiag])
d1         = ufl.as_vector([-offdiag, diag])

m = u[0] - 0.5 - kappa1/ufl.pi*ufl.atan(kappa2*u[1])
s = ufl.as_vector([dt/tau * u[0] *  (1.-u[0]) * m , u[0]])

a = ( alpha*alpha*dt/tau *
        ( ufl.inner( ufl.dot(d0, ufl.grad(u[0])), ufl.grad(v[0])[0]) +
          ufl.inner( ufl.dot(d1, ufl.grad(u[0])), ufl.grad(v[0])[1]) ) +
      2.25*dt * ufl.inner(ufl.grad(u[1]), ufl.grad(v[1])) +
      ufl.inner(u,v) -
      ufl.inner(s,v)
    ) * ufl.dx(0)

model.generate(a)
# now add implicit part of forcing to source and un coefficient function
lhsModel = model.makeAndImport(grid2d,name="crystal_left").get()


# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
initial_gf  =  grid2d.globalGridFunction("initial", initial)
solution    = sp.interpolate(initial_gf, name="solution")
solution_n  = sp.interpolate(initial_gf, name="solution_n")
forcing     = sp.interpolate([0,0,0], name="forcing")
# left hand side scheme
solver    = scheme.create("FemScheme", solution, lhsModel, "left")
# right hand side scheme
rhs       = scheme.create("FemScheme", forcing,  rhsModel, "rhs")

lhsModel.setun(solution_n)

# start adaptation
hgrid = grid2d.hierarchicalGrid

marker = hgrid.marker
def mark(element):
    # return marker.keep
    solutionLocal = solution.localFunction(element)
    x = [1/3,1/3]
    grad = solutionLocal.jacobian(x)
    if grad[0].infinity_norm > 1.:
      return marker.refine if element.level < maxLevel+1 else marker.keep
    else:
      return marker.coarsen

for i in range(0,maxLevel+1):
    hgrid.mark(mark)
    hgrid.adapt([solution])
    hgrid.loadBalance([solution])
    solution.interpolate(initial_gf)

solution_n.assign(solution)

# time lopp
# ---------
count    = 0
t        = 0.
savestep = saveinterval
grid2d.writeVTK("crystal", pointdata=[solution], celldata=[level_gf], number=count)

while t < endTime:
    rhs(solution_n, forcing)
    solver.solve(target=solution, rhs=forcing)
    t += dt
    print('count: ',count,"t = ",t)
    if t > savestep:
        savestep += saveinterval
        count += 1
        grid2d.writeVTK("crystal", pointdata=[solution], celldata=[level_gf], number=count)
    hgrid.mark(mark)
    hgrid.adapt([solution])
    hgrid.loadBalance([solution])
    solution_n.assign(solution)

print("END")
