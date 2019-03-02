import math
from ufl import *

from dune.grid import cartesianDomain
from dune.ufl import Space

import dune.create as create

grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
spc = create.space("Lagrange", grid, dimrange=1, order=2, storage="fem")

uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

d = 0.001
p = 1.7

rhs = (x[0] + x[1]) * v[0]

a = (pow(d + inner(grad(u), grad(u)), (p-2)/2)*inner(grad(u), grad(v)) + inner(u, v)) * dx + 10*inner(u, v) * ds
#b = sin(2*math.pi*x[0])*sin(2*math.pi*x[1]) * v[0] * dx
b = rhs * dx + 10*rhs * ds

model = create.model("integrands", grid, a==b)

scheme = create.scheme("galerkin", model, spc)

uh = spc.interpolate([0],name="solution")
scheme.solve(target=uh)
grid.writeVTK("p-laplace", pointdata=[uh])
