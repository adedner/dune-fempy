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
grid2d    = fem.leafGrid("../data/sphere.dgf", "ALUSimplexGrid", dimgrid=2, dimworld=3)
# discrete function for Gamma(t) and setup surface grid
print("positions")
positions = grid2d.interpolate(lambda x: x, space="Lagrange", name="positions", variant="global")
surface   = gridpart.create("Geometry", positions )
# vtk writer
vtk       = surface.vtkWriter()
# space for discrete solution on Gamma(t)
sp        = space.create( "Lagrange", surface, dimrange=3, polorder=1 )
# final time and time step
endTime = 1.
dt      = 0.001


# set up left and right hand side models
# --------------------------------------
model     = duneuflmodel.DuneUFLModel(3,3)
u         = model.trialFunction()
v         = model.testFunction()

a = ( dt*0.5*ufl.inner(ufl.grad(u),ufl.grad(v)) +
      ufl.inner(u,v) ) * ufl.dx(0)
model.generate(a,"mcf_left")
mcfModel = model.makeAndImport(surface).get()

a = ( -dt*0.5*ufl.inner(ufl.grad(u),ufl.grad(v)) +
      ufl.inner(u,v) ) * ufl.dx(0)
model.generate(a,"mcf_right")
rhsModel = model.makeAndImport(surface).get()

# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
solution  = sp.interpolate( lambda x: x, name="solution", variant="global" )
forcing   = sp.interpolate( [0,0,0], name="forcing" )
# left hand side scheme
solver    = scheme.create( "FemScheme", solution, mcfModel, "mcf" )
# right hand side scheme
rhs       = scheme.create( "FemScheme", forcing,  rhsModel, "rhs" )
solution.addToVTKWriter(vtk, vtk.PointData)
forcing.addToVTKWriter(vtk, vtk.PointData)

# time lopp
# ---------
count   = 0
t       = 0.
vtk.write("mcf"+str(count));

while t<endTime:
    forcing.clear()
    rhs(solution,forcing)
    solution.clear()
    solver.solve(target=solution, rhs=forcing)
    t     += dt
    count += 1
    vtk.write("mcf"+str(count))
    positions.assign( solution.dofVector() )
