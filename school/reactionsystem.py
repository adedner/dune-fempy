import math
import ufl
from ufl import grad, div, dot, dx, inner, sin, cos, pi
import dune.ufl
import dune.grid
import dune.fem

endTime  = 0.2
saveInterval = 0.1

gridView = dune.grid.structuredGrid([0,0],[1,1],[100,100])

# Define stream-function and transport field w
# analytic:
# Psi     = 0.1*sin(2*pi*x[0])*sin(2*pi*x[1])
# discrete:
velocitySpace = dune.fem.space.lagrange(gridView, order=1)
Psi  = velocitySpace.interpolate(0,name="stream_function")
u    = ufl.TrialFunction(velocitySpace)
phi  = ufl.TestFunction(velocitySpace)
x    = ufl.SpatialCoordinate(velocitySpace)
form = ( inner(grad(u),grad(phi)) -
         0.1*2*(2*pi)**2*sin(2*pi*x[0])*sin(2*pi*x[1]) * phi ) * dx
dbc  = dune.ufl.DirichletBC(velocitySpace,0)
velocityScheme = dune.fem.scheme.galerkin([form == 0, dbc], solver="cg")
velocityScheme.solve(target=Psi)
w       = ufl.as_vector([-Psi.dx(1),Psi.dx(0)])

space = dune.fem.space.lagrange(gridView, order=1, dimRange=3)
u     = ufl.TrialFunction(space)
phi   = ufl.TestFunction(space)
x     = ufl.SpatialCoordinate(space)

# reaction, diffusion and other coefficients
K   = dune.ufl.Constant(10.0, "reactionRate")
eps = dune.ufl.Constant(0.01, "diffusionRate")
dt = dune.ufl.Constant(0.01, "timeStep")
t  = dune.ufl.Constant(0., "timeStep")

# define storage for discrete solutions
uh     = space.interpolate([0,0,0], name="uh")
uh_old = uh.copy()

# define source terms
Q  = dune.ufl.Constant(0.1, "sourceStrength")
P1 = ufl.as_vector([0.1,0.1])
P2 = ufl.as_vector([0.9,0.9])
RF = 0.075
f1 = ufl.conditional(dot(x-P1,x-P1) < RF**2, Q, 0)
f2 = ufl.conditional(dot(x-P2,x-P2) < RF**2, Q, 0)
f  = ufl.as_vector([f1,f2,0])

# reaction rates
r   = K*ufl.as_vector([u[0]*u[1], u[0]*u[1], -2*u[0]*u[1] + 10*u[2]])

xForm = (dot(grad(u)*w + r - f, phi) + eps * inner(grad(u), grad(phi))) * dx
form  = dot(u - uh_old, phi) * dx + dt * xForm

scheme = dune.fem.scheme.galerkin(form == 0, solver="gmres")

nextSaveTime = saveInterval
vtk = gridView.sequencedVTK("reactionsystem", pointdata=[uh], pointvector={"velocity":w})
vtk()

while t.value < endTime:
    uh_old.assign(uh)
    info = scheme.solve(target=uh)
    t.value += dt.value
    print("Computed solution at time",t.value,\
              "iterations: ", info["linear_iterations"], "#Ent: ", gridView.size(0) )
    if t.value >= nextSaveTime or t.value >= endTime:
        vtk()
        nextSaveTime += saveInterval
