"""Solve the Laplace equation
"""

from mpi4py import MPI

import math
from ufl import *

from dune.models.elliptic import importModel
import dune.ufl
import dune.fem
import dune.fem.function as gf

# dune.fem.create.spaceGenerator.force = True

dune.femmpi.parameter.append("../data/parameter")

dune.femmpi.parameter.append("../data/parameter")

grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
spc = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1, field="double")
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
imag = Coefficient(uflSpace)
const = VectorConstant(triangle,2)
coeff = Coefficient(uflSpace)
const0 = Constant(triangle)
x = SpatialCoordinate(uflSpace.cell())

f = const0*const[0]*coeff[0]* (cos(2*math.pi*x[0])*cos(2*math.pi*x[1]) + cos(2*math.pi*x[0])*cos(2*math.pi*x[1]))

a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
b = f * v[0] * dx

model = importModel(grid, a==b).get()
model.setConstant(const, [10.,-10.])
model.setConstant(const0, [20.]) # note: still using a list here instead of a double
gfunc = gf.MathExpression(["1."])
coeffFunc = grid.globalGridFunction("global_velocity", gfunc)
model.setCoefficient(coeff, coeffFunc)

scheme = dune.fem.create.scheme("FemScheme", spc, model,\
        "scheme", storage="Istl",\
   parameters={"fem.solver.newton.linabstol": 1e-9,
    "fem.solver.newton.linreduction": 1e-9,
    "fem.solver.newton.verbose": 1,
    "fem.solver.newton.linear.verbose": 1})

grid.writeVTK("laplace", pointdata=[ scheme.solve() ])
