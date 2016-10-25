from __future__ import print_function

import math
from ufl import *

from dune.ufl import Space as UFLSpace
from dune.source import SourceWriter
from dune.models.integrands import compileUFL

uflSpace = UFLSpace(2, 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

a = inner(grad(u), grad(v)) * dx
a += u[0]*inner(u, v) * dx
a += inner(avg(grad(u)), jump(grad(v))) * dS

b = v[0] * ds

integrands = compileUFL(a == b, tempVars=True)

# write model to file
# -------------------

writer = SourceWriter("myintegrands.hh")
writer.openNameSpace('demo')
integrands.write(writer, "MyModel")
writer.closeNameSpace('demo')
