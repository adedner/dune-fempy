from __future__ import print_function

import math
from ufl import *

from dune.grid import cartesianDomain
from dune.alugrid import aluConformGrid

from dune.ufl import Space as UFLSpace
from dune.source import SourceWriter
from dune.source.cplusplus import NameSpace
from dune.models.integrands import compileUFL, load

import dune.create as create

domain = cartesianDomain([0, 0], [1, 1], [8, 8])
grid = aluConformGrid(domain, dimgrid=2)
space = create.space("Lagrange", grid, dimrange=1, order=2, storage="istl")

uflSpace = UFLSpace(2, 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

a = inner(u, v) * dx
b = sin(math.pi*x[0]) * sin(math.pi * x[1]) * v[0] * dx

integrands = compileUFL(a == b)


# write model to file
# -------------------

code = NameSpace('demo')
code.append(integrands.code("MyIntegrands"))

SourceWriter("myintegrands.hh").emit(code)


# load integrands
# ---------------

module = load(grid, integrands)
integrands = module.Integrands()

scheme = create.scheme("galerkin", space, integrands)
solution, _ = scheme.solve()
grid.writeVTK("l2projection", pointdata=[solution])
