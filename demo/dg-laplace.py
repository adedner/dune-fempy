from __future__ import print_function

import math
from ufl import *

from dune.grid import cartesianDomain
from dune.fem import parameter
from dune.ufl import Space

import dune.create as create

parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
spc = create.space("dgonb", grid, dimrange=1, order=2, storage="istl")

uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
n = FacetNormal(uflSpace.cell())
mu = 7.5 * 16

a = inner(grad(u), grad(v)) * dx
a -= (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
a += mu * inner(jump(u), jump(v)) * dS
a -= (inner(outer(u, n), grad(v)) + inner(grad(u), outer(v, n))) * ds
a += mu * inner(u, v) * ds

b = sin(pi*x[0])*sin(pi*x[1])*v[0]*dx

model = create.model("integrands", grid, a == b)

newtonParameter = {"linabstol": 1e-13, "linreduction": 1e-13, "tolerance": 1e-12, "verbose": "true", "linear.verbose": "false"}
scheme = create.scheme("galerkin", spc, model, parameters={"fem.solver.newton." + k: v for k, v in newtonParameter.items()})

solution, _ = scheme.solve()
grid.writeVTK("dg-laplace", celldata=[solution])
