from __future__ import print_function

from mpi4py import MPI

import math
import ufl

import dune.ufl
import dune.models.elliptic

import dune.fem

# Basic setup
# -----------
# set up reference domain
grid = dune.fem.leafGrid("../data/sphere.dgf", "ALUSimplexGrid", dimgrid=2, dimworld=3)
grid.hierarchicalGrid.globalRefine(2)

# discrete function for Gamma(t) and setup surface grid
positions = grid.interpolate(lambda x: x * (1.+0.5*math.sin(2.*math.pi*x[0]*x[1])*math.cos(math.pi*x[2])),
            space="Lagrange", name="positions")
surface   = dune.fem.create.gridpart("Geometry", positions)
# space for discrete solution on Gamma(t)
spc       = dune.fem.create.space("Lagrange", surface, dimrange=surface.dimWorld, polorder=1)
# final time and time step
endTime   = 0.25
dt        = 0.0025
theta     = 0.5

# set up left and right hand side models
# --------------------------------------
print("surface.dimWorld: ", surface.dimWorld)
uflSpace = dune.ufl.Space(surface.dimGrid, surface.dimWorld, surface.dimWorld)
# uflSpace = dune.ufl.Space(3,3)
u   = ufl.TrialFunction(uflSpace)
v   = ufl.TestFunction(uflSpace)
u_n = ufl.Coefficient(uflSpace)

a_im = (dt * theta * ufl.inner(ufl.grad(u), ufl.grad(v)) + ufl.inner(u, v)) * ufl.dx(0)
a_ex = (-dt * (1-theta) * ufl.inner(ufl.grad(u), ufl.grad(v)) + ufl.inner(u, v)) * ufl.dx(0)
lhsModel = dune.models.elliptic.importModel(surface, dune.models.elliptic.compileUFL(a_im == 0)).get()
rhsModel = dune.models.elliptic.importModel(surface, dune.models.elliptic.compileUFL(a_ex == 0)).get()

# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
solution   = spc.interpolate(lambda x: x, name="solution")
forcing = spc.interpolate([0,]*surface.dimWorld, name="solution")
# left hand side scheme
solver    = dune.fem.create.scheme("FemScheme", solution, lhsModel, "lhs")
rhs       = dune.fem.create.scheme("FemScheme", solution, rhsModel, "rhs")

# time loop
# ---------
count   = 0
t       = 0.
surface.writeVTK("mcf", pointdata=[solution], number=count)

while t < endTime:
    rhs(solution, forcing)
    solver.solve(forcing, solution )
    t     += dt
    count += 1
    surface.writeVTK("mcf", pointdata=[solution], number=count)
    positions.assign(solution.dofVector())
