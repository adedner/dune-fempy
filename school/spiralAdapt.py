# <markdowncell>
# # Spiral Wave
#
# This demonstrates the simulation of spiral waves in an excitable media. It consists of system of reaction diffusion equations with two components. Both the model parameters and the approach for discretizing the system are taken from http://www.scholarpedia.org/article/Barkley_model.
#
# We use the _Barkley model_ in its simplest form:
# \begin{align*}
#   \frac{\partial u}{\partial_t}
#        &= \frac{1}{\varepsilon}f(u,v) + \Delta u \\
#   \frac{\partial v}{\partial_t} &= h(u,v)
# \end{align*}
# where
# \begin{gather}
#   f(u,v(=u\Big(1-u\Big)\Big(u-\frac{v+b}{a}\Big)
# \end{gather}
# The function $h$ can take different forms, e.g., in its simplest form
# \begin{gather}
#   h(u,v) = u - v~.
# \end{gather}
# Finally, $\varepsilon,a,b$ for more details on how to chose these parameters check the web page provided above.
#
# We employ a carefully constructed linear time stepping scheme for this model: let $u^n,v^n$ be given functions approximating the solution at a time $t^n$. To compute approximations $u^{m+1},v^{m+1}$ at a later time
# $t^{n+1}=t^n+\tau$ we first split up the non linear function $f$ as follows:
# \begin{align*}
#   f(u,v) = f_I(u,u,v) + f_E(u,v)
# \end{align*}
# where using $u^*(V):=\frac{V+b}{a}$:
# \begin{align*}
#   f_I(u,U,V) &= \begin{cases}
#     u\;(1-U)\;(\;U-U^*(V)\;) & U < U^*(V) \\
#     -u\;U\;(\;U-U^*(V)\;)    & U \geq U^*(V)
#   \end{cases} \\
# \text{and} \\
#     f_E(U,V) &= \begin{cases}
#     0 & U < U^*(V) \\
#     U\;(\;U-U^*(V)\;)    & U \geq U^*(V)
#   \end{cases} \\
# \end{align*}
# Thus $f_I(u,U,V) = -m(U,V)u$ with
# \begin{align*}
#   m(U,V) &= \begin{cases}
#     (U-1)\;(\;U-U^*(V)\;) & U < U^*(V) \\
#     U\;(\;U-U^*(V)\;)    & U \geq U^*(V)
#   \end{cases}
# \end{align*}
# Note that $u,v$ are assumed to take values only between zero and one so that therefore $m(u^n,v^n) > 0$. Therefore, the following time discrete version of the Barkley model has a linear, positive definite elliptic operator on its left hand side:
# \begin{align*}
#   -\tau\Delta u^{n+1} +
#    (1+\frac{\tau}{\varepsilon} m(u^n,v^n))\; u^{n+1}
#        &= u^n + \frac{\tau}{\varepsilon} f_E(u^n,v^n) \\
#   v^{n+1} &= v^n + \tau h(u^n,v^n)
# \end{align*}
# Which can now be solved using a finite element discretization for $u^n,v^n$.
#
# Note that by taking the slow reaction $h(u,v)$ explicitly, the equation for $v^{n+1}$ is purely algebraic. We will therefore construct a scalar model for computing $u^{n+1}$ only and compute $v^{{n+1}}$ be using the interpolation method on the space applied to
# $v^n + \tau h(u^n,v^n)$.
#
# Let's get started by importing some standard python packages, ufl, and some part of the dune-fempy package:
# <codecell>

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
maxLevel     = 13
startLevel   = 9
dt           = dune.ufl.Constant(0.5,"dt")
t            = dune.ufl.Constant(0,"time")
endTime      = 7.
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

domain   = dune.grid.cartesianDomain([0,0],[2.5,2.5],[5,5])
baseView = dune.alugrid.aluConformGrid(domain)
gridView = dune.fem.view.adaptiveLeafGridView( baseView )
gridView.hierarchicalGrid.globalRefine(startLevel)

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

solverParameters =\
       {"newton.tolerance": 1e-10,
        "newton.linear.tolerance": 1e-8,
        "newton.verbose": False,
        "newton.linear.verbose": False}
scheme = dune.fem.scheme.galerkin( equation, solver="cg", parameters=solverParameters)

# <markdowncell>
# Error estimator
# <codecell>

fvspace = dune.fem.space.finiteVolume(uh.space.grid)
estimate = fvspace.interpolate([0], name="estimate")

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
        gridView.writeVTK("spiral", pointdata=[uh,vh], number=count, celldata=[estimate,levelFunction])
        nextSaveTime += saveInterval
        count += 1
    if t.value > 5:
        dune.fem.mark(estimate, maxTol, 0.1 * maxTol, 0,maxLevel)
        dune.fem.adapt([uh,vh])
