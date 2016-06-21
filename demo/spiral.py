from __future__ import print_function
import math
from mpi4py import MPI

import ufl
import dune.models.femufl as duneuflmodel
import dune.fem as fem
import dune.fem.gridpart as gridpart
import dune.fem.space as space
import dune.fem.scheme as scheme

# http://www.scholarpedia.org/article/Barkley_model
endTime    = 20.
dt         = 0.05
spiral_a   = 0.75
spiral_b   = 0.02
spiral_eps = 0.02
dimRange   = 2

def initial(x):
    return [ 1   if x[1]>1.75 else 0,\
             0.5 if x[0]<1.75 else 0. ]

#################################################################
### Extra model code:
implicitCode = """\
      double a=0.75;
      double b=0.02;
      double eps=0.02;
      double dt=0.025;
      RangeType unValue;
      unLocal_->evaluate( point, unValue );
      double uth = (unValue[1]+b)/a;
      if ( unValue[0] <= uth )
        flux[ 0 ] += -dt/eps * value[0] * (1.-unValue[0]) * (unValue[0]-uth);
      else
        flux[ 0 ] +=  dt/eps * value[0] * unValue[0] * (unValue[0]-uth);
"""
explicitCode = """\
      double a=0.75;
      double b=0.02;
      double eps=0.02;
      double dt=0.025;
      double uth = (value[1]+b)/a;
      if (value[0]>uth)
        flux[0] += dt/eps * value[0] * (value[0]-uth);
"""
#################################################################

# Basic setup
# -----------
# set up reference domain
grid2d    = fem.leafGrid("../data/spiral-2d.dgf", "YaspGrid", dimgrid=2, dimworld=2)
# vtk writer
vtk       = grid2d.vtkWriter()
sp        = space.create( "Lagrange", grid2d, dimrange=dimRange, polorder=1 )

# set up left and right hand side models
# --------------------------------------
model     = duneuflmodel.DuneUFLModel(2,dimRange)
u         = model.trialFunction()
v         = model.testFunction()

# right hand sie (time derivative part + explicit forcing in v)
a = ( ufl.inner(u,v) + dt*ufl.inner(u[0]-u[1],v[1]) ) * ufl.dx(0)
model.generate(a,"spiral_right" )
# now add implicit part of forcing to source
model.add2Source( explicitCode )
rhsModel = model.makeAndImport(grid2d).get()

# left hand side (heat equation in first variable + backward Euler in time)
a = ( dt/100.*ufl.inner(ufl.grad(u[0]),ufl.grad(v[0])) +
      ufl.inner(u,v) ) * ufl.dx(0)
model.generate(a,"spiral_left")
# now add implicit part of forcing to source and un coefficient function
model.addCoefficient( "un", 2 )
model.add2Source( implicitCode )
lhsModel = model.makeAndImport(grid2d).get()


# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
solution    = sp.interpolate( initial, name="solution", variant="global" )
solution_n  = sp.interpolate( initial, name="solution_n", variant="global" )
forcing     = sp.interpolate( [0,0,0], name="forcing" )
# left hand side scheme
solver    = scheme.create( "FemScheme", solution, lhsModel, "left" )
# right hand side scheme
rhs       = scheme.create( "FemScheme", forcing,  rhsModel, "rhs" )
solution.addToVTKWriter(vtk, vtk.PointData)


lhsModel.setun(solution_n)

# time lopp
# ---------
count   = 0
t       = 0.
vtk.write("spiral"+str(count));

while t<endTime:
    rhs(solution,forcing)
    solver.solve(target=solution, rhs=forcing)
    t     += dt
    count += 1
    vtk.write("spiral"+str(count))
    solution_n.assign(solution)

print("END")
