from __future__ import print_function
import math
from mpi4py import MPI

import ufl
import dune.models.femufl as duneuflmodel
import dune.fem as fem
import dune.fem.gridpart as gridpart
import dune.fem.space as space
import dune.fem.scheme as scheme

# Basic setup
# -----------
# set up reference domain
grid2d = fem.leafGrid("../data/sphere.dgf", "ALUSimplexGrid", dimgrid=2, dimworld=3)
grid2d.hierarchicalGrid.globalRefine(2)

# discrete function for Gamma(t) and setup surface grid
positions = grid2d.interpolate(lambda x:
                x * (1.+0.5*math.sin(2.*math.pi*x[0]*x[1])*math.cos(math.pi*x[2])),
            space="Lagrange", name="positions")
surface   = gridpart.create("Geometry", positions)
# space for discrete solution on Gamma(t)
sp        = space.create( "Lagrange", surface, dimrange=3, polorder=1)
# final time and time step
endTime   = 0.25
dt        = 0.0025


# set up left and right hand side models
# --------------------------------------
model     = duneuflmodel.DuneUFLModel(3,3)
u         = model.trialFunction()
v         = model.testFunction()

a = ( dt*0.5*ufl.inner(ufl.grad(u),ufl.grad(v)) +
      ufl.inner(u,v) ) * ufl.dx(0)
model.generate(a)
mcfModel = model.makeAndImport(surface,name="mcf_left").get()

model.clear()

a = ( -dt*0.5*ufl.inner(ufl.grad(u),ufl.grad(v)) +
      ufl.inner(u,v) ) * ufl.dx(0)
model.generate(a)
rhsModel = model.makeAndImport(surface,name="mcf_right").get()

# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
solution  = sp.interpolate(lambda x: x, name="solution")
forcing   = sp.interpolate([0,0,0], name="forcing" )
# left hand side scheme
solver    = scheme.create("FemScheme", solution, mcfModel, "mcf")
# right hand side scheme
rhs       = scheme.create("FemScheme", forcing,  rhsModel, "rhs")

# time lopp
# ---------
count   = 0
t       = 0.
surface.writeVTK("mcf", pointdata=[solution], number=count)

while t<endTime:
    rhs(solution, forcing)
    solver.solve(forcing, solution)
    t     += dt
    count += 1
    surface.writeVTK("mcf", pointdata=[solution], number=count)
    positions.assign(solution.dofVector())
