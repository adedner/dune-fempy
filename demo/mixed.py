import dune.fem
from dune.fem.plotting import plotPointData as plot
import dune.create as create
from dune.ufl import DirichletBC, Space

## simple example of mixed BCs (u[0] with zero Neumann, u[1] with Dirichlet = 2)
dune.fem.parameter.append("parameter")

dimDomain = 2
dimRange = 2
grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [4, 4]), dimgrid=2)

from ufl import TestFunction, TrialFunction, SpatialCoordinate
uflSpace = dune.ufl.Space(dimDomain, dimRange)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

from math import pi,log10
from ufl import cos, as_vector, dx, grad, inner
f = (8*pi*pi+1)*cos(2*pi*x[0])*cos(2*pi*x[1])
equation = (inner(grad(u), grad(v)) + inner(u,v)) * dx == f * v[0] * dx

spc = create.space("Lagrange", grid, dimrange=dimRange, order=1)

boundary_cds = [None, 2]
model = create.model("elliptic", grid, equation, DirichletBC(uflSpace, boundary_cds, 6))
scheme = create.scheme("h1", spc, model)
solution, _ = scheme.solve()
