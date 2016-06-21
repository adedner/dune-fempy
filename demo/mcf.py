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
positions = grid2d.interpolate(lambda x: x, space="Lagrange", name="psoitions", variant="global")
surface   = gridpart.create("Geometry", positions )
# vtk writer
vtk       = surface.vtkWriter()
# space for discrete solution on Gamma(t)
sp        = space.create( "Lagrange", surface, dimrange=3, polorder=1 )
# final time and time step
endTime = 1.
dt      = 0.01


# set up left and right hand side models
# --------------------------------------
model     = duneuflmodel.DuneUFLModel(2,3)
u         = model.trialFunction()
v         = model.testFunction()
x         = model.spatialCoordinate()
tau       = model.coefficient('tau')
old       = model.coefficient('old')

model.setCoefficient("tau",[dt])   # set time step
a = ( tau[0]*ufl.inner(ufl.grad(u),ufl.grad(v)) + ufl.inner(u,v) ) * ufl.dx(0)
model.generate(a,"mcf_left")
Model = model.makeAndImport(surface)
mcfModel  = Model.get()

a = ( ufl.inner(u,v)) * ufl.dx(0)
model.generate(a,"mcf_right")
Model = model.makeAndImport(surface)
rhsModel  = Model.get()

# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
solution  = sp.interpolate( lambda x: x, variant="global" )
forcing   = sp.interpolate( lambda x: x, variant="global" )
# left hand side scheme
solver    = scheme.create( "FemScheme", solution, mcfModel, "mcf" )
# right hand side scheme
rhs       = scheme.create( "FemScheme", forcing, rhsModel, "rhs" )

# time lopp
# ---------
count   = 0
t       = 0.
solution.addToVTKWriter(vtk, vtk.PointData)
vtk.write("mcf"+str(count));

while t<endTime:
    rhs(solution,forcing)
    solver.solve(target=solution, rhs=forcing)
    t     += dt
    count += 1
    vtk.write("mcf"+str(count))
    positions.assign( solution.dofVector() )
