from __future__ import print_function, division

import math
import numpy

from ufl import *

from dune.grid import cartesianDomain
from dune.fem import parameter
from dune.ufl import Space

import dune.create as create

parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu-0", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

vertices = numpy.array([(0,0), (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1)])
triangles = numpy.array([(0,1,2), (0,2,3), (0,3,4), (0,4,5), (0,5,6), (0,6,7)])

grid = create.grid("ALUConform", {"vertex": vertices, "simplex": triangles}, dimgrid=2)
grid.hierarchicalGrid.globalRefine(4)

spc = create.space("Lagrange", grid, dimrange=1, order=1, storage="istl")

uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*math.pi, 0)
exact = as_vector([inner(x,x)**(0.5*180/270) * sin((180/270) * phi)])

model = create.model("elliptic", grid, inner(grad(u), grad(v))*dx == 0, dirichlet={1: exact})

newtonParameter = {"linabstol": 1e-13, "linreduction": 1e-13, "tolerance": 1e-12, "verbose": "true", "linear.verbose": "false"}
scheme = create.scheme("h1", spc, model, parameters={"fem.solver.newton." + k: v for k, v in newtonParameter.items()})

solution, _ = scheme.solve()

fvspc = create.space("FiniteVolume", grid, dimrange=1, storage="istl")
estimate = fvspc.interpolate(solution, name="estimate")

hT = MaxCellEdgeLength(uflSpace.cell())
he = MaxFacetEdgeLength(uflSpace.cell())
n = FacetNormal(uflSpace.cell())
estimator_ufl = hT**2 * (div(grad(u[0])))**2 * v[0] * dx + he * inner(jump(grad(u[0])), n('+'))**2 * avg(v[0]) * dS
estimator_model = create.model("integrands", grid, estimator_ufl == 0)
estimator = create.operator("galerkin", estimator_model, spc, fvspc)

estimator(solution, estimate)

exact_grid = create.function("ufl", grid, "exact", 2, exact)
grid.writeVTK("reentrant-corner", pointdata=[solution, exact_grid], celldata=[estimate])
