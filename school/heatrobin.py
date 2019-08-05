import ufl
from ufl import grad, div, dot, dx, ds, inner, sin, cos, pi, exp, sqrt
import dune.ufl
import dune.grid
import dune.fem

endTime  = 0.1
saveInterval = 1 # for VTK

gridView = dune.grid.structuredGrid([-1,-1],[1,1],[40,40])

space = dune.fem.space.lagrange(gridView, order=1)
u     = ufl.TrialFunction(space)
phi   = ufl.TestFunction(space)
x     = ufl.SpatialCoordinate(space)
dt    = dune.ufl.Constant(5e-2, "timeStep")
t     = dune.ufl.Constant(0.0, "time")

# define storage for discrete solutions
uh     = space.interpolate(0, name="uh")
uh_old = uh.copy()

# initial solution
initial = 0

# problem definition

# moving oven
ROven = 0.6
omegaOven = 0.01 * pi * t
P = ufl.as_vector([ROven*cos(omegaOven*t), ROven*sin(omegaOven*t)])
rOven = 0.2
ovenEnergy = 8
chiOven = ufl.conditional(dot(x-P, x-P) < rOven**2, 1, 0)
ovenLoad = ovenEnergy * chiOven

# desk in corner of room
deskCenter = [-0.8, -0.8]
deskSize = 0.2
chiDesk = ufl.conditional(abs(x[0]-deskCenter[0]) < deskSize, 1, 0)\
  * ufl.conditional(abs(x[1] - deskCenter[1]) < deskSize, 1, 0)

# Robin condition for window
windowWidth = 0.5
transmissionCoefficient = 1.2
outerTemperature = -5.0
chiWindow = ufl.conditional(abs(x[1]-1.0) < 1e-8, 1, 0)*ufl.conditional(abs(x[0]) < windowWidth, 1, 0)
rBC = transmissionCoefficient * (u - outerTemperature) * chiWindow

# heat diffussion
K = 0.01

# space form
diffusiveFlux = K*grad(u)
source = -ovenLoad
xForm = (dot(diffusiveFlux, grad(phi)) + source * phi) * dx + rBC * phi * ds

# add time discretization
form = dot(u - uh_old, phi) * dx + dt * xForm

scheme = dune.fem.scheme.galerkin(form == 0, solver="cg")

nextSaveTime = saveInterval

uh.interpolate(initial)

vtk = gridView.sequencedVTK("heatrobin", pointdata=[uh])
vtk()

while t.value < endTime:
    uh_old.assign(uh)
    info = scheme.solve(target=uh)
    t.value += dt.value
    deskTemperature = dune.fem.function.integrate(gridView, uh * chiDesk, order=1) / deskSize**2 / 4

    print("Computed solution at time", t.value,
              "desk temperature", deskTemperature,
              "iterations: ", info["linear_iterations"],
              "#Ent: ", gridView.size(0) )
    if t.value >= nextSaveTime or t.value >= endTime:
        vtk()
        nextSaveTime += saveInterval
