from __future__ import print_function
from mpi4py import MPI

import math
from ufl import *
import dune.ufl
import dune.models.elliptic

uflSpace = dune.ufl.Space(2, 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

a = (inner(grad(u), grad(v)) + inner(u, v)) * dx(0)
b = sin(2*math.pi*x[0])*sin(2*math.pi*x[1]) * v[0] * dx(0) - x[0]*(1-x[0]) * v[0] * ds(0)

model = dune.models.elliptic.compileUFL(a == b, dimRange=1)

writer = dune.models.elliptic.SourceWriter("mymodel.hh")
writer.openNameSpace('demo')
model.write(writer, "MyModel")
writer.closeNameSpace('demo')
