from __future__ import print_function
import math
from mpi4py import MPI

import ufl
import dune.models.femufl as duneuflmodel
import dune.fem as fem
import dune.fem.gridpart as gridpart
import dune.fem.space as space
import dune.fem.scheme as scheme

#################################################################
dimRange     = 2
gamma        = 0.675
eps          = 0.015
dt           = 0.001
endTime      = 10
saveinterval = 0.02
def initial(x):
    r  = (x-[6,6]).two_norm
    r0 = 0.5;
    return [ (math.tanh(-(r-r0)/0.1)+1.)*0.5 , -0.5 ]

#################################################################

# Basic setup
# -----------
# set up reference domain
grid2d    = fem.leafGrid("../data/crystal-2d.dgf", "ALUSimplexGrid", dimgrid=2, dimworld=2, refinement="conforming")
grid2d.hierarchicalGrid.globalRefine(1)
# vtk writer
vtk       = grid2d.vtkWriter()
sp        = space.create( "Lagrange", grid2d, dimrange=dimRange, polorder=1 )

# set up left and right hand side models
# --------------------------------------
model     = duneuflmodel.DuneUFLModel(2,dimRange)
u         = model.trialFunction()
v         = model.testFunction()

# right hand sie (time derivative part + explicit forcing in v)
a = ( ufl.inner(u,v) - ufl.inner(u[0],v[1]) ) * ufl.dx(0)
model.generate(a, name="crystal_right" )
# now add implicit part of forcing to source
rhsModel = model.makeAndImport(grid2d).get()
model.clear()

# left hand side (heat equation in first variable + backward Euler in time)
dun = model.coefficient('dun',2)
N          = 6.
psi        = ufl.pi/8. + ufl.atan( dun[1] / (dun[0]+1e-8) )
Phi        = ufl.tan(N/2.*psi)
beta       = ( 1 - Phi*Phi )/ (1 + Phi*Phi)
dbeta_dPhi = -2. * N * Phi / (1 + Phi*Phi)
fac        = 1 + 0.02*beta
diag       = fac*fac
offdiag    = -fac * 0.02 * dbeta_dPhi
d0         = ufl.as_vector( [ diag, offdiag ] )
d1         = ufl.as_vector( [ -offdiag, diag ] )

c = 0.5 + 0.9/ufl.pi*ufl.atan(20.*u[1])
s = ufl.as_vector( [ u[0] *  (1.-u[0]) * (u[0]-c) * gamma/eps/eps*dt , u[0] ] )

a = ( gamma*dt *
        ( ufl.inner( ufl.dot( d0, ufl.grad(u[0]) ), ufl.grad(v[0])[0] ) +
          ufl.inner( ufl.dot( d1, ufl.grad(u[0]) ), ufl.grad(v[0])[1] ) ) +
      2.25*dt * ufl.inner( ufl.grad(u[1]), ufl.grad(v[1]) ) +
      ufl.inner(u,v) -
      ufl.inner(s,v)
    ) * ufl.dx(0)

model.generate(a, name="crystal_left")
# now add implicit part of forcing to source and un coefficient function
lhsModel = model.makeAndImport(grid2d).get()


# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
initial_gf  =  grid2d.globalGridFunction("initial", initial)
solution    = sp.interpolate( initial_gf, name="solution" )
solution_n  = sp.interpolate( initial_gf, name="solution_n" )
forcing     = sp.interpolate( [0,0,0], name="forcing" )
# left hand side scheme
solver    = scheme.create( "FemScheme", solution, lhsModel, "left" )
# right hand side scheme
rhs       = scheme.create( "FemScheme", forcing,  rhsModel, "rhs" )
solution.addToVTKWriter(vtk, vtk.PointData)

class LocalExprA:
    def __init__(self,df):
        self.df = df
        self.dimR = 2
def dunLocal(en,x):
  jac = solution_n.localFunction(en).jacobian(x)
  return [ jac[0][0],jac[0][1] ]
dun_gf = grid2d.localGridFunction( "nabla_un0", dunLocal )

lhsModel.setdun(dun_gf)

# start adaptation
maxLevel = 6
hgrid = grid2d.hierarchicalGrid

marker = hgrid.marker
def mark(element):
    # return marker.keep
    solutionLocal = solution.localFunction(element)
    x = [1/3,1/3]
    grad = solutionLocal.jacobian( x )
    if grad[0].infinity_norm > 0.1:
      return marker.refine if element.level < maxLevel+1 else marker.keep
    else:
      return marker.coarsen

for i in range(0,maxLevel+1):
    hgrid.mark(mark)
    hgrid.adapt( [solution] )
    hgrid.loadBalance( [solution] )
    solution.interpolate( initial_gf )

solution_n.assign( solution )

# time lopp
# ---------
count    = 0
t        = 0.
savestep = saveinterval
vtk.write( "crystal"+str(count) )

while t < endTime:
    rhs( solution_n, forcing )
    solver.solve( target=solution, rhs=forcing )
    t     += dt
    print('count: ',count,"t = ",t)
    if t > savestep:
        vtk.write( "crystal"+str(count) )
        count += 1
        savestep += saveinterval

    hgrid.mark( mark )
    hgrid.adapt( [solution] )
    hgrid.loadBalance( [solution] )
    solution_n.assign( solution )

print("END")
