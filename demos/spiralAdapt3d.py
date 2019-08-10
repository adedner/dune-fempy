import math
import ufl
from ufl import div, grad, inner, dot, dx, dS, jump, avg
import dune.ufl
import dune.grid
import dune.fem
import dune.alugrid

from dune.grid import cartesianDomain, gridFunction

# <markdowncell>
# In our attempt we will discretize the model as a 2x2 system. Here are some possible model parameters and initial conditions (we even have two sets of model parameters to choose from):
# <codecell>

linearSpiral = True
maxLevel     = 13
startLevel   = 3
dt           = dune.ufl.Constant(0.1,"dt")
t            = dune.ufl.Constant(0,"time")
endTime      = 20.
saveInterval = 0.5
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
# Now we set up the reference domain, the Lagrange finite element space (second order), and discrete functions for $(u^n,v^n($, $(u^{n+1},v^{n+1})$:
# <codecell>

#domain   = dune.grid.cartesianDomain([0,0],[2.5,2.5],[5,5])
domain   = dune.grid.cartesianDomain([0,0,0],[2.5,2.5,2.5],[5,5,5])
#baseView = dune.alugrid.aluConformGrid(domain)
baseView = dune.grid.ugGrid(domain)
gridView = dune.fem.view.adaptiveLeafGridView( baseView )
gridView.hierarchicalGrid.globalRefine(startLevel)

maxOrder = 1
space = dune.fem.space.lagrange( gridView, order=maxOrder,   storage="istl" )
spcpm = dune.fem.space.lagrange( gridView, order=maxOrder-1, storage="istl" )

x = ufl.SpatialCoordinate(space)
iu = lambda s: ufl.conditional(s > 1.25, 1, 0 )
top = ufl.conditional( x[2] > 1.25,1,0)
initial_u = iu(x[1])*top + iu(2.5-x[1])*(1.0 - top)
#initial_u = ufl.conditional( x[2] > 1.25, iu(x[1]), iu(2.5 - x[1]))
initial_v = ufl.conditional(x[0]<1.25,0.5,0)

uh     = space.interpolate( initial_u, name="u" )
uh_pm1 = spcpm.interpolate( initial_u, name="u_p-1" )
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

def markp(element):
    return 2

# initial mark
#dune.fem.spaceAdapt(space,markp,[uh])

def markpm1(element):
    return space.localOrder( element ) - 1


# <markdowncell>
# The model is now completely implemented and can be created, together with the corresponding scheme:
# <codecell>

solverParameters =\
       {"newton.tolerance": 1e-10,
        "newton.linear.tolerance": 1e-8,
        "newton.linear.preconditioning.method": "amg-ilu",
        "newton.verbose": False,
        "newton.linear.verbose": False}
scheme = dune.fem.scheme.galerkin( equation, solver="cg", parameters=solverParameters)

# <markdowncell>
# Error estimator
# <codecell>

fvspace = dune.fem.space.finiteVolume(uh.space.grid)
estimate = fvspace.interpolate([0], name="estimate")
estimate_pm1 = fvspace.interpolate([0], name="estimate_pm1")

chi = ufl.TestFunction(fvspace)
hT  = ufl.MaxCellEdgeLength(fvspace.cell())
he  = ufl.MaxFacetEdgeLength(fvspace.cell())('+')
n   = ufl.FacetNormal(fvspace.cell())

residual = (u-uh_n)/dt - div(diffusiveFlux) + source(u,u,u,vh)

estimator_ufl = hT**2 * residual**2 * chi * dx +\
                he * inner( jump(diffusiveFlux), n('+'))**2 * avg(chi) * dS
estimator = dune.fem.operator.galerkin(estimator_ufl)

# <markdowncell>
# Time loop
# <codecell>

nextSaveTime = saveInterval
count = 0
levelFunction = dune.fem.function.levelFunction(gridView)
gridView.writeVTK("spiral", pointdata=[uh,vh], number=count, celldata=[estimate,levelFunction])
count += 1

@gridFunction(gridView,name="pEstimate")
def pEstimator(e,x):
    r    = estimate.localFunction(e).evaluate(x)
    r_p1 = estimate_pm1.localFunction(e).evaluate(x)
    # if r[0] < 1e-15:
    #     eta = 0
    # else:
    #     eta = abs(r[0]-r_p1[0]) / r[0]
    eta = abs(sum(r)-sum(r_p1))
    return [eta, sum(r), sum(r_p1)]

def markpDiff(element):
    # ptol = 1e-5    # for relative
    ptol = 1e-16 # 1e-17   # absolute
    # evaluate smoothness indicator
    eta = pEstimator(element,element.referenceElement.center)
    # get current polorder
    polorder = space.localOrder(element)

    newPolorder = polorder
    if eta[0] < ptol:
        newPolorder = polorder-1 if polorder > 1 else polorder
    elif eta[0] > 100.*ptol:
        newPolorder = polorder+1 if polorder < maxOrder else polorder

    return newPolorder

@gridFunction(gridView,"pDegree")
def pDegree(element,x):
    return space.localOrder(element)

while t.value < endTime:
    uh_n.assign(uh)
    vh_n.assign(vh)
    info = scheme.solve(target=uh)
    t.value += dt.value
    print("Computed solution at time", t.value,
          "iterations: ", info["linear_iterations"],
          "#Ent: ", gridView.size(0) )
    vh.interpolate( vh_n + dt*spiral_h(uh_n, vh_n) )

    estimator(uh, estimate)
    maxEst = max(estimate.dofVector)
    print("max est: ", maxEst)
    if t.value >= nextSaveTime or t.value >= endTime:
        gridView.writeVTK("spiral", pointdata=[uh,vh], number=count, celldata=[estimate,levelFunction,pDegree])
        nextSaveTime += saveInterval
        count += 1
    if t.value > 5:
        dune.fem.mark(estimate, maxTol, 0.1 * maxTol, 0,maxLevel)
        dune.fem.adapt([uh,vh])

        #estimator(uh, estimate)
        # mark space with one p lower locally
        #dune.fem.spaceAdapt(spcpm,markpm1,[uh_pm1])
        #uh_pm1.interpolate( uh )
        #estimator(uh_pm1, estimate_pm1)
        #dune.fem.spaceAdapt(space,markpDiff, [uh])
