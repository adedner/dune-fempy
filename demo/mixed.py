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
from ufl import cos, sin,as_vector, dx, grad, inner
f = (8*pi*pi+1)*cos(2*pi*x[0])*cos(2*pi*x[1])
exact = as_vector( [(8*pi*pi+1)*cos(2*pi*x[0])*cos(2*pi*x[1]),x[0]*x[1] ] )
equation = (inner(grad(u), grad(v)) + inner(u,v)) * dx == f * v[0] * dx

spc = create.space("Lagrange", grid, dimrange=dimRange, order=1)

model = create.model("elliptic", grid, equation,
        DirichletBC(uflSpace, [None,2], 2),
        DirichletBC(uflSpace, [x[0]*x[1],None], 4),
        DirichletBC(uflSpace, [sin(x[0]),cos(x[1])], 6),
        DirichletBC(uflSpace, exact, 8))
scheme = create.scheme("h1", spc, model)
solution, _ = scheme.solve()
