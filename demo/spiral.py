from __future__ import print_function
import math
from mpi4py import MPI
from functools import reduce

import ufl
import dune.models.femufl as duneuflmodel
import dune.fem as fem
import dune.fem.gridpart as gridpart
import dune.fem.space as space
import dune.fem.scheme as scheme

# http://www.scholarpedia.org/article/Barkley_model
dimRange   = 2
endTime    = 30.
dt         = 0.1
if 1:
    spiral_a   = 0.75
    spiral_b   = 0.02
    spiral_eps = 0.02
    spiral_D   = 1./100
    def spiral_h(u,v): return u - v
else:
    spiral_a   = 0.75
    spiral_b   = 0.0006
    spiral_eps = 0.08
    def spiral_h(u,v): return u**3 - v

def initial(x):
    return [ 1   if x[1]>1.75 else 0,\
             0.5 if x[0]<1.75 else 0. ]

#################################################################
### Extra model code:
sourceCode = """\
      RangeType unValue;
      unLocal_->evaluate( point, unValue );
      double uth = (unValue[1]+spiral_b)/spiral_a;
      if ( unValue[0] <= uth )
        flux[ 0 ] -= dt/spiral_eps * value[0] * (1.-unValue[0]) * (unValue[0]-uth);
      else
        flux[ 0 ] -= dt/spiral_eps * unValue[0] * (1.-value[0]) * (unValue[0]-uth);
"""
linSourceCode = """\
      RangeType unValue;
      unLocal_->evaluate( point, unValue );
      double uth = (unValue[1]+spiral_b)/spiral_a;
      if ( unValue[0] <= uth )
        flux[ 0 ] -= dt/spiral_eps * value[0] * (1.-unValue[0]) * (unValue[0]-uth);
      else
        flux[ 0 ] -= dt/spiral_eps * unValue[0] * (-value[0]) * (unValue[0]-uth);
"""
repls = ('spiral_a', str(spiral_a)), ('spiral_b', str(spiral_b)),\
        ('spiral_eps', str(spiral_eps)), ('dt',str(dt))
sourceCode    = reduce(lambda a, kv: a.replace(*kv), repls, sourceCode)
linSourceCode = reduce(lambda a, kv: a.replace(*kv), repls, linSourceCode)
#################################################################

# Basic setup
# -----------
# set up reference domain
grid2d    = fem.leafGrid("../data/spiral-2d.dgf", "YaspGrid", dimgrid=2, dimworld=2)
sp        = space.create( "Lagrange", grid2d, dimrange=dimRange, polorder=1 )

# set up left and right hand side models
# --------------------------------------
ufl2model = duneuflmodel.DuneUFLModel(grid2d.dimWorld,dimRange)
u         = ufl2model.trialFunction()
v         = ufl2model.testFunction()
un        = ufl2model.coefficient('un',dimRange)

# right hand sie (time derivative part + explicit forcing in v)
a_ex = ( ufl.inner(un,v) + dt*ufl.inner( spiral_h( un[0],un[1]), v[1] ) ) * ufl.dx(0)
# left hand side (heat equation in first variable + backward Euler in time)
a_im = ( dt*spiral_D*ufl.inner(ufl.grad(u[0]),ufl.grad(v[0])) +
         ufl.inner(u,v) ) * ufl.dx(0)
ufl2model.generate(a_im-a_ex)
# now add implicit part of forcing to source and un coefficient function
ufl2model.add2Source( sourceCode, linSourceCode )
model = ufl2model.makeAndImport(grid2d,name="spiral").get()

# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
solution    = sp.interpolate( initial, name="solution" )
solution_n  = sp.interpolate( initial, name="solution_n" )
forcing     = sp.interpolate( [0,0,0], name="forcing" )
solver    = scheme.create( "FemScheme", solution, model, "spiral" )

model.setun(solution_n)

# time lopp
# ---------
count   = 0
t       = 0.
grid2d.writeVTK("spiral", pointdata=[solution], number=count)

while t<endTime:
    solver.solve( target=solution, assemble=(count==0) )
    t     += dt
    count += 1
    grid2d.writeVTK("spiral", pointdata=[solution], number=count)
    solution_n.assign( solution )

print("END")
