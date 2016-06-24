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
## www.ctcms.nist.gov/fipy/examples/phase/generated/examples.phase.anisotropy.html
dimRange     = 2
dt           = 5.e-4
endTime      = 1
saveinterval = 0.001
maxLevel     = 8
def initial(x):
    r  = (x-[6,6]).two_norm
    r0 = 0.1
    return [ (math.tanh(-(r-r0)/0.1)+1.)*0.5 , -0.5 ]

#################################################################

# Basic setup
# -----------
# set up reference domain
grid2d    = fem.leafGrid("../data/crystal-2d.dgf", "ALUSimplexGrid", dimgrid=2, dimworld=2, refinement="conforming")
grid2d.hierarchicalGrid.globalRefine(3)
sp        = space.create( "Lagrange", grid2d, dimrange=dimRange, polorder=1 )
level_gf  = grid2d.localGridFunction("level", function.Levels())

# set up left and right hand side models
# --------------------------------------
model      = duneuflmodel.DuneUFLModel(2,dimRange)
rhsModel   = model.makeAndImport(grid2d,name="crystal_right",header="crystal_rightModel.hh").get()
model.addCoefficient("dun", dimRange)
lhsModel   = model.makeAndImport(grid2d,name="crystal_left",header="crystal_leftModel.hh").get()


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
grid2d.writeVTK( "crystal", pointdata=[solution], celldata=[level_gf], number=count )

while t < endTime:
    rhs( solution_n, forcing )
    solver.solve( target=solution, rhs=forcing )
    t     += dt
    print('count: ',count,"t = ",t)
    if t > savestep:
        savestep += saveinterval
        count += 1
        grid2d.writeVTK( "crystal", pointdata=[solution], celldata=[level_gf], number=count )
    hgrid.mark( mark )
    hgrid.adapt( [solution] )
    hgrid.loadBalance( [solution] )
    solution_n.assign( solution )

print("END")
