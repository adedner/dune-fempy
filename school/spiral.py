import math
import ufl
from ufl import div, grad, inner, dot, dx, dS, jump, avg
import dune.ufl
import dune.grid
import dune.fem
import dune.alugrid

# <markdowncell>
# In our attempt we will discretize the model as a 2x2 system. Here are some possible model parameters and initial conditions (we even have two sets of model parameters to choose from):
# <codecell>

linearSpiral = True
dt           = dune.ufl.Constant(0.1,"dt")
t            = dune.ufl.Constant(0,"time")
endTime      = 1
saveInterval = 0.5
maxTol       = 1e-6

if linearSpiral:
    spiral_a   = 0.75
    spiral_b   = 0.02
    spiral_eps = 0.02
    spiral_D   = 1./100
    def spiral_h(u,v): return u - v
else:
    spiral_a   = 0.75
    spiral_b   = 0.0006
    spiral_eps = 0.08
    def spiral_h(u,v): return u**3 - v


# <markdowncell>
# Now we set up the reference domain, the Lagrange finite element space (second order), and discrete functions for $(u^n,v^n($, $(u^{n+1},v^{n+1})$:
# <codecell>

domain   = dune.grid.cartesianDomain([0,0],[2.5,2.5],[100,100])
gridView = dune.alugrid.aluConformGrid(domain)

space = dune.fem.space.lagrange( gridView, order=1 )

x = ufl.SpatialCoordinate(space)
initial_u = ufl.conditional(x[1]>1.25,1,0)
initial_v = ufl.conditional(x[0]<1.25,0.5,0)

uh   = space.interpolate( initial_u, name="u" )
uh_n = uh.copy()
vh   = space.interpolate( initial_v, name="v" )
vh_n = vh.copy()


# <markdowncell>
# Setting up the model
# <codecell>

u   = ufl.TrialFunction(space)
phi = ufl.TestFunction(space)

ustar          = lambda v: (v+spiral_b)/spiral_a

diffusiveFlux  = spiral_D * grad(u)
source         = lambda u1,u2,u3,v: -1/spiral_eps * u1*(1-u2)*(u3-ustar(v))

xForm  = inner(diffusiveFlux, grad(phi)) * dx
xForm += ufl.conditional(uh_n<ustar(vh_n), source(u,uh_n,uh_n,vh_n), source(uh_n,u,uh_n,vh_n)) * phi * dx
form   = ( inner(u,phi) - inner(uh_n, phi) ) * dx + dt*xForm

equation   = form == 0

# <markdowncell>
# The model is now completely implemented and can be created, together with the corresponding scheme:
# <codecell>

scheme = dune.fem.scheme.galerkin( equation, solver="cg" )

residual = (u-uh_n)/dt - div(diffusiveFlux) + source(u,u,u,vh)

# <markdowncell>
# Time loop
# <codecell>

nextSaveTime = saveInterval

vtk = gridView.sequencedVTK("spiral", pointdata=[uh,vh])
vtk()

while t.value < endTime:
    uh_n.assign(uh)
    vh_n.assign(vh)
    info = scheme.solve(target=uh)
    vh.interpolate( vh_n + dt*spiral_h(uh_n, vh_n) )
    t.value += dt.value
    print("Computed solution at time", t.value,
          "iterations: ", info["linear_iterations"],
          "#Ent: ", gridView.size(0) )

    if t.value >= nextSaveTime or t.value >= endTime:
        vtk()
        nextSaveTime += saveInterval
