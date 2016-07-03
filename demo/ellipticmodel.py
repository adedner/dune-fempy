from __future__ import print_function

import math
from ufl import *

from dune.ufl import Space as UFLSpace
from dune.models.elliptic import compileUFL, SourceWriter

uflSpace = UFLSpace(2, 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

a = (inner(grad(u), grad(v)) + u[0]*u[0]*v[0]) * dx(0)
b = sin(2*math.pi*x[0])*sin(2*math.pi*x[1]) * v[0] * dx(0) - x[0]*(1-x[0]) * v[0] * ds(0)

model = compileUFL(a == b, dirichlet={1:[x[0]], 2:[x[1]], 3:[zero(tuple())]}, tempVars = False)

# write model to file
# -------------------

writer = SourceWriter("mymodel.hh")
writer.openNameSpace('demo')
model.write(writer, "MyModel")
writer.closeNameSpace('demo')
