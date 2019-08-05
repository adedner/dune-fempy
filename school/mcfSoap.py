import ufl
from ufl import grad, dot, dx, inner, Identity, sqrt
import dune.ufl
import dune.grid
import dune.fem
from dune.alugrid import aluConformGrid as leafGridView
from dune.fem.view import geometryGridView
from dune.fem.space import lagrange as solutionSpace
from dune.fem.scheme import galerkin as solutionScheme

# setup reference surface
referenceView = leafGridView("soap.dgf", dimgrid=2, dimworld=3)
referenceView.hierarchicalGrid.globalRefine(8)
space = solutionSpace(referenceView, dimRange=referenceView.dimWorld, order=1)

# setup deformed surface
x         = ufl.SpatialCoordinate(space)
positions = space.interpolate(x/sqrt(dot(x,x)), name="position")
gridView  = geometryGridView(positions)
space     = solutionSpace(gridView, dimRange=gridView.dimWorld, order=1)

u     = ufl.TrialFunction(space)
phi   = ufl.TestFunction(space)
dt    = dune.ufl.Constant(0.05, "timeStep")
t     = dune.ufl.Constant(0.0, "time")

# define storage for discrete solutions
uh     = space.interpolate(x, name="uh")
uh_old = uh.copy()

# problem definition

# space form
xForm = inner(grad(u), grad(phi)) * dx

# add time discretization
form = dot(u - uh_old, phi) * dx + dt * xForm

# define dirichlet boundary conditions
bc = dune.ufl.DirichletBC(space,x)

# setup scheme
scheme = solutionScheme([form == 0,bc], space, solver="cg")

endTime      = 0.1
saveInterval = 0.05
nextSaveTime = saveInterval

vtk = gridView.sequencedVTK("mcfSoap", pointdata=[uh])
vtk()

while t.value < endTime:
    uh_old.assign(uh)
    info = scheme.solve(target=uh)
    t.value += dt.value

    positions.dofVector.assign(uh.dofVector)

    print("Computed solution at time", t.value,
              "iterations: ", info["linear_iterations"],
              "#Ent: ", gridView.size(0) )
    if t.value >= nextSaveTime or t.value >= endTime:
        vtk()
        nextSaveTime += saveInterval
