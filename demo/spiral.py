from __future__ import print_function

import math
import ufl

import dune.ufl
import dune.models.elliptic
import dune.fem

import dune.create as create

from functools import reduce

# http://www.scholarpedia.org/article/Barkley_model
dimRange   = 2
endTime    = 3. # 30.
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

# Basic setup
# -----------
# set up reference domain
grid = create.grid("Yasp", "../data/spiral-2d.dgf", dimgrid=2)
spc  = create.space( "Lagrange", grid, dimrange=dimRange, order=1 )

# set up left and right hand side models
# --------------------------------------
uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), dimRange)
u = ufl.TrialFunction(uflSpace)
v = ufl.TestFunction(uflSpace)
un = ufl.Coefficient(uflSpace)

# right hand sie (time derivative part + explicit forcing in v)
a_ex = (ufl.inner(un, v) + dt * ufl.inner(spiral_h(un[0], un[1]), v[1])) * ufl.dx
# left hand side (heat equation in first variable + backward Euler in time)
a_im = (dt * spiral_D * ufl.inner(ufl.grad(u[0]), ufl.grad(v[0])) + ufl.inner(u,v)) * ufl.dx

modelCode = dune.models.elliptic.compileUFL(a_im == a_ex)

# extra model source code
# -----------------------
sourceCode = """\
RangeType un;
coefficient< UNCOUNT >().evaluate( x, un );
double uth = (un[ 1 ] + spiral_b) / spiral_a;
if( un[ 0 ] <= uth )
  result[ 0 ] -= dt/spiral_eps * u[ 0 ] * (1.0 - un[ 0 ]) * (un[ 0 ] - uth);
else
  result[ 0 ] -= dt/spiral_eps * un[ 0 ] * (1.0 - u[ 0 ]) * (un[ 0 ] - uth);
"""
linSourceCode = """\
RangeType un;
coefficient< UNCOUNT >().evaluate( x, un );
double uth = (un[ 1 ] + spiral_b) / spiral_a;
if( un[ 0 ] <= uth )
  result[ 0 ] -= dt/spiral_eps * u[ 0 ] * (1.0 - un[ 0 ]) * (un[ 0 ] - uth);
else
  result[ 0 ] -= dt/spiral_eps * un[ 0 ] * (-u[ 0 ]) * (un[ 0 ] - uth);
"""
repls = ('spiral_a', str(spiral_a)), ('spiral_b', str(spiral_b)),\
        ('spiral_eps', str(spiral_eps)), ('dt',str(dt)),\
        ('UNCOUNT', "0")

modelCode.source.append(reduce(lambda a, kv: a.replace(*kv), repls, sourceCode))
modelCode.linSource.append(reduce(lambda a, kv: a.replace(*kv), repls, linSourceCode))

# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
solution    = spc.interpolate( initial, name="solution" )
solution_n  = spc.interpolate( initial, name="solution_n" )
forcing     = spc.interpolate( [0,0,0], name="forcing" )

model = create.model("elliptic", grid, modelCode, coefficients={un:solution_n} )

scheme = create.scheme("h1", solution, model)

# time loop
# ---------
count   = 0
t       = 0.
grid.writeVTK("spiral", pointdata=[solution], number=count)

while t < endTime:
    print(">>> Computing solution a t = " + str(t + dt))
    solution_n.assign(solution)
    scheme.solve(target=solution) # , assemble=(count==0))
    t     += dt
    count += 1
    grid.writeVTK("spiral", pointdata=[solution], number=count)

grid = None
print("END")
