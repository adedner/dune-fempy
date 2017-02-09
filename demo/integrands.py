from __future__ import print_function

import math
from ufl import *

from dune.ufl import Space as UFLSpace
from dune.source import SourceWriter
from dune.source.cplusplus import NameSpace
from dune.models.integrands import compileUFL

uflSpace = UFLSpace(2, 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
n = FacetNormal(uflSpace.cell())

print(outer(u, n))

a = inner(grad(u), grad(v)) * dx
a += u[0]*inner(u, v) * dx
a -= inner(outer(u('+'), n('-')), grad(v)('-')) * dS

b = v[0] * ds

integrands, _ = compileUFL(a == b, tempVars=True)


# write model to file
# -------------------

code = NameSpace('demo')
code.append(integrands.code("MyIntegrands"))

SourceWriter("myintegrands.hh").emit(code)


# load integrands
# ---------------

from dune.grid import cartesianDomain
from dune.alugrid import aluConformGrid
from dune.models.integrands import create

domain = cartesianDomain([0, 0], [1, 1], [8, 8])
grid = aluConformGrid(domain, dimgrid=2)

integrands = create(grid, integrands)
