from __future__ import print_function

from mpi4py import MPI

import math
from ufl import *

from dune.ufl import Space as UFLSpace
from dune.models.elliptic import compileUFL, SourceWriter

uflSpace = UFLSpace(2, 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

a = inner(grad(u),grad(v)) * dx
H = grad(grad(u[0]))
a = a + H[0,0] * v[0] * dx
a = a + det(H) * v[0] * dx

b = v[0] * ds

model = compileUFL(a == b, dirichlet={1:[x[0]], 2:[x[1]], 3:[zero(tuple())]}, tempVars = False)

# write model to file
# -------------------

writer = SourceWriter("mymodel.hh")
writer.openNameSpace('demo')
model.write(writer, "MyModel")
writer.closeNameSpace('demo')
