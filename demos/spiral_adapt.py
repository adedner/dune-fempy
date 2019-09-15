# <markdowncell>
# In this example we revisit the Barkley model for spiral waves and as before
# we have two sets of model parameters to choose from. We want to add
# dynamic grid adaptivity for different grid structures based on a residual
# based error estimator similar to the one used for elliptic problems in
# the previous example. To allow for both conforming and non conforming
# refinement (i.e. to handle hanging nodes) we use a DG type stabilization
# approach. Note that since we will be using a Lagrange space these
# skeleton terms evaluate to zero over conforming intersections.
# <codecell>

import sys
import math
import ufl
from ufl import div, grad, inner, dot, dx, dS, jump, avg, outer
import dune.ufl
import dune.grid
import dune.fem
import dune.alugrid

from dune.grid import cartesianDomain, gridFunction
from dune.fem import parameter

# <markdowncell>
# We can choose between a 2d and a 3d setup and use both conforming and
# non-conforming refinement (in the first case a bisection simplex grid is
# used in the second case a simplex grid with quartering/subdivision into
# eight children will be used. We can also fix the order (default is 2).
# As in the non adaptive example we have two possible choices for the right
# hand side function $h$ to choose from.
# <codecell>

dim = 3
if len(sys.argv) > 1:
    dim = int(sys.argv[1])
print("Using dim = ",dim)
nonConforming = True

order = 2
if len(sys.argv) > 2:
    maxOrder = int(sys.argv[2])
print("Using order = ",order)

linearSpiral = True

# <markdowncell>
# Next we setup the constant used in the model and for the discretization:
# <codecell>

maxLevel     = 5 if dim == 2 else 4
startLevel   = 2 if dim == 2 else 1
dt           = dune.ufl.Constant(0.1,"dt")
t            = dune.ufl.Constant(0,"time")
endTime      = 15.
saveInterval = 1.0
maxTol       = 1e-4

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
# Now we set up the reference domain, the Lagrange finite element space, and discrete functions for $(u^n,v^n($, $(u^{n+1},v^{n+1})$:
# Note: for the adaptive simulations we need to use the `adaptiveLeafGridView`
# <codecell>

domain   = dune.grid.cartesianDomain([0,]*dim,[2.5,]*dim,[5,]*dim)
if nonConforming:
    baseView = dune.alugrid.aluSimplexGrid(domain)
else:
    baseView = dune.alugrid.aluConformGrid(domain)

gridView = dune.fem.view.adaptiveLeafGridView( baseView )

maxLevel *= gridView.hierarchicalGrid.refineStepsForHalf
startLevel *= gridView.hierarchicalGrid.refineStepsForHalf
gridView.hierarchicalGrid.globalRefine(startLevel)

print("Initial grid size = ",gridView.size( 0 ))

space = dune.fem.space.lagrangehp( gridView, order=order, storage="istl" )

x = ufl.SpatialCoordinate(space)
iu = lambda s: ufl.conditional(s > 1.25, 1, 0 )
top = 1.0
if dim == 3:
    top = ufl.conditional( x[2] > 1.25,1,0)

initial_u = iu(x[1])*top + iu(2.5-x[1])*(1.0 - top)
initial_v = ufl.conditional(x[0]<1.25,0.5,0)

uh   = space.interpolate( initial_u, name="u" )
uh_n = uh.copy()
vh   = space.interpolate( initial_v, name="v" )
vh_n = vh.copy()

u   = ufl.TrialFunction(space)
phi = ufl.TestFunction(space)
hT  = ufl.MaxCellEdgeLength(space.cell())
hS  = ufl.avg( ufl.MaxFacetEdgeLength(space.cell()) )
hs =  ufl.MaxFacetEdgeLength(space.cell())('+')
n   = ufl.FacetNormal(space.cell())

ustar          = lambda v: (v+spiral_b)/spiral_a
diffusiveFlux  = lambda w,d: spiral_D * d
source         = lambda u1,u2,u3,v: -1/spiral_eps * u1*(1-u2)*(u3-ustar(v))
source         = lambda u1,u2,u3,v: -1/spiral_eps * u1*(1-u2)*(u3-ustar(v))

xForm  = inner(diffusiveFlux(u,grad(u)), grad(phi)) * dx
xForm += ufl.conditional(uh_n<ustar(vh_n), source(u,uh_n,uh_n,vh_n), source(uh_n,u,uh_n,vh_n)) * phi * dx

# <markdowncell>
# To handle the possible non conforming intersection we add DG type
# skeleton terms:
# <codecell>

penalty = 5 * (order * ( order + 1 )) * spiral_D
if nonConforming:
    xForm -= ( inner( outer(jump(u), n('+')), avg(diffusiveFlux(u,grad(phi)))) +\
               inner( avg(diffusiveFlux(u,grad(u))), outer(jump(phi), n('+'))) ) * dS
    xForm += penalty/hS * inner(jump(u), jump(phi)) * dS

# <markdowncell>
# After adding the time discretization (a simple backward Euler scheme),
# the model is now completely implemented and the scheme can be created:
# <codecell>

form   = ( inner(u,phi) - inner(uh_n, phi) ) * dx + dt*xForm
solverParameters =\
       {"newton.tolerance": 1e-8,
        "newton.linear.tolerance": 1e-12,
        "newton.linear.preconditioning.method": "amg-ilu",
        "newton.linear.maxiterations":1000,
        "newton.verbose": False,
        "newton.linear.verbose": True}
scheme = dune.fem.scheme.galerkin( form==0, solver="cg", parameters=solverParameters)

# <markdowncell>
# Next we define the error estimator based on a residual type indicator
# similar to the one used for the Laplace problem:
# <codecell>

fvspace = dune.fem.space.finiteVolume(uh.space.grid)
estimate = fvspace.interpolate([0], name="estimate")

chi = ufl.TestFunction(fvspace)

residual = (u-uh_n)/dt - div(diffusiveFlux(u,grad(u))) + source(u,u,u,vh)

estimator_ufl = hT**2 * residual**2 * chi * dx +\
                hS * inner( jump(diffusiveFlux(u,grad(u))), n('+'))**2 * avg(chi) * dS +\
                1/hS * jump(u)**2 * avg(chi) * dS
estimator = dune.fem.operator.galerkin(estimator_ufl)

# <markdowncell>
# Now the time loop - note that we wait until the spiral has developed from
# the initial conditions before the adaptive refinement/coarsening kicks
# off:
# <codecell>
nextSaveTime = saveInterval
levelFunction = dune.fem.function.levelFunction(gridView)
subSampling = 1 if order < 3 else 2
vtk = gridView.sequencedVTK("spiral", pointdata=[uh,vh], subsampling=subSampling)
vtk_grid = gridView.sequencedVTK("spiral-grid", pointdata=[uh,vh], celldata=[estimate,levelFunction])
vtk()
vtk_grid()

while t.value < endTime:
    uh_n.assign(uh)
    vh_n.assign(vh)
    info = scheme.solve(target=uh)
    t.value += dt.value
    vh.interpolate( vh_n + dt*spiral_h(uh_n, vh_n) )

    estimator(uh, estimate)
    maxEst = max(estimate.dofVector)
    print("max est: ", maxEst)
    if t.value >= nextSaveTime-0.01 or t.value >= endTime:
        print("Writing vtu at time ", t.value,
              "#Ent: ", gridView.size(0) )
        vtk()
        vtk_grid()
        nextSaveTime = t.value + saveInterval
    if t.value > 5: # start the adaptation once the spiral has developed
        dune.fem.mark(estimate, maxTol, 0.05 * maxTol, 0,maxLevel)
        dune.fem.adapt([uh,vh])
