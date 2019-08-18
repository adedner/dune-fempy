import ufl
from ufl import grad, dot, dx, inner, Identity, sqrt, sin, pi, cos
import dune.ufl
import dune.grid
import dune.fem
from dune.alugrid import aluConformGrid as leafGridView
from dune.fem.view import geometryGridView
from dune.fem.space import lagrange as solutionSpace
from dune.fem.scheme import galerkin as solutionScheme

order = 3
storage = "istl"
# setup reference surface
referenceView = leafGridView("sphere.dgf", dimgrid=2, dimworld=3)
space = solutionSpace(referenceView, dimRange=referenceView.dimWorld,
        order=order, storage=storage)

# setup deformed surface
x         = ufl.SpatialCoordinate(space)
# positions = space.interpolate(x, name="position")
positions = space.interpolate(x * (1 + 0.5*sin(2*pi*(x[0]+x[1]))*cos(0.25*pi*x[2])), name="position")
gridView  = geometryGridView(positions)
space     = solutionSpace(gridView, dimRange=gridView.dimWorld,
            order=order, storage=storage)

u     = ufl.TrialFunction(space)
phi   = ufl.TestFunction(space)
dt    = dune.ufl.Constant(0.01, "timeStep")
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
dune.fem.parameter.append({"fem.verboserank": 0})
solverParameters =\
       {"newton.tolerance": 1e-9,
        "newton.linear.tolerance": 1e-11,
        "newton.linear.preconditioning.method": "ilu",
        "newton.verbose": False,
        "newton.linear.verbose": False}
scheme = solutionScheme([form == 0], space, solver="cg",
                        parameters=solverParameters)

endTime      = 0.3
saveInterval = 0.01
nextSaveTime = saveInterval

vtk = gridView.sequencedVTK("mcfPaper", pointdata=[uh], subsampling=2)
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
