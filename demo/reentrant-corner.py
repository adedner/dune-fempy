from __future__ import print_function, division

import math
import numpy

from ufl import *

from dune.grid import cartesianDomain, Marker
from dune.fem import parameter, adapt
from dune.fem.function import levelFunction
from dune.ufl import Space

import dune.create as create

parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

vertices = numpy.array([(0,0), (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1)])
triangles = numpy.array([(2,1,0), (0,3,2), (4,3,0,), (0,5,4), (6,5,0), (0,7,6)])
domain = {"vertices": vertices, "simplices": triangles}
grid   = create.view("adaptive", grid="ALUConform", constructor=domain, dimgrid=2)
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

fvspc = create.space("finitevolume", grid, dimrange=1, storage="istl")
estimate = fvspc.interpolate([0], name="estimate")

hT = MaxCellEdgeLength(uflSpace.cell())
he = MaxFacetEdgeLength(uflSpace.cell())('+')
n = FacetNormal(uflSpace.cell())
estimator_ufl = hT**2 * (div(grad(u[0])))**2 * v[0] * dx + he * inner(jump(grad(u[0])), n('+'))**2 * avg(v[0]) * dS
estimator_model = create.model("integrands", grid, estimator_ufl == 0)
estimator = create.operator("galerkin", estimator_model, spc, fvspc)

tolerance = 0.001
gridSize = grid.size(0)
def mark(element):
    estLocal = estimate(element, element.geometry.referenceElement.center)
    return Marker.refine if estLocal[0] > tolerance / gridSize else Marker.keep

for i in range(20):
    solution, _ = scheme.solve()
    estimator(solution, estimate)
    grid.writeVTK("reentrant-corner",
            pointdata={"solution":solution, "exact":exact},
            celldata=[estimate,levelFunction(grid)], number=i)
    eta = sum(estimate.dofVector)
    marked = grid.hierarchicalGrid.mark(mark)
    print(gridSize, eta, marked)
    if eta < tolerance or sum(marked)==0:
        break
    adapt(grid.hierarchicalGrid,[solution])
    gridSize = grid.size(0)
