import math
import ufl
from ufl import grad, dot, dx, inner, Identity, sqrt, sin, cos, pi
import dune.ufl
import dune.grid
import dune.fem
from dune.alugrid import aluConformGrid as leafGridView
from dune.fem.view import geometryGridView
from dune.fem.space import lagrange as solutionSpace
from dune.fem.scheme import galerkin as solutionScheme

endTime      = 0.05
saveInterval = 0.05

# setup reference surface
referenceView = leafGridView("sphere.dgf", dimgrid=2, dimworld=3)
space = solutionSpace(referenceView, dimRange=referenceView.dimWorld, order=1)

# setup deformed surface
x         = ufl.SpatialCoordinate(space)
# positions = space.interpolate(x, name="position")
positions = space.interpolate(x*(1 + 0.5*sin(2*pi*x[0]*x[1])*cos(pi*x[2])), name="position")
gridView  = geometryGridView(positions)
space     = solutionSpace(gridView, dimRange=gridView.dimWorld, order=1)

u     = ufl.TrialFunction(space)
phi   = ufl.TestFunction(space)
dt    = dune.ufl.Constant(0.001, "timeStep")
t     = dune.ufl.Constant(0.0, "time")

# define storage for discrete solutions
uh     = space.interpolate(x, name="uh")
uh_old = uh.copy()

# problem definition

# space form
xForm = inner(grad(u), grad(phi)) * dx
# add time discretization
form = dot(u - uh_old, phi) * dx + dt * xForm

# setup scheme
scheme = solutionScheme(form == 0, space, solver="cg")

nextSaveTime = saveInterval

vtk = gridView.sequencedVTK("mcfSphere", pointdata=[uh])
vtk()

while t.value < endTime:
    uh_old.assign(uh)
    info = scheme.solve(target=uh)
    t.value += dt.value

    positions.dofVector.assign(uh.dofVector)

    [radius,area] = dune.fem.function.integrate(gridView,[sqrt(dot(x,x)),1],order=4)

    print("Computed solution at time", t.value,
              "error in radius: ",radius/area- math.sqrt(1-4*t.value),
              "iterations: ", info["linear_iterations"],
              "#Ent: ", gridView.size(0) )
    if t.value >= nextSaveTime or t.value >= endTime:
        vtk()
        nextSaveTime += saveInterval
