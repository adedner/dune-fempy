from __future__ import print_function

from mpi4py import MPI

import math
from ufl import *

import dune.models.elliptic
import dune.ufl
import dune.grid
import dune.fem
import dune.fem.space
import dune.fem.scheme

grid = dune.grid.create("ALUConform", dune.grid.cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
spc = dune.fem.space.create("Lagrange", grid, dimrange=1, order=2)

uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

d = 0.001
p = 1.7

rhs = (x[0] + x[1]) * v[0]

a = (pow(d + inner(grad(u), grad(u)), (p-2)/2)*inner(grad(u), grad(v)) + inner(u, v)) * dx + 10*inner(u, v) * ds
#b = sin(2*math.pi*x[0])*sin(2*math.pi*x[1]) * v[0] * dx
b = rhs * dx + 10*rhs * ds

Model = dune.fem.create.ellipticModel(grid,a==b)

scheme = dune.fem.scheme.create("h1", spc, Model(), "scheme")
grid.writeVTK("p-laplace", pointdata=[scheme.solve()])
