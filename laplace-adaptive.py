# <markdowncell>
# ## Re-entrant Corner Problem
#
# Here we will consider the classic _re-entrant corner_ problem,
# \begin{align*}
#   -\Delta u &= f, && \text{in } \Omega, \\
#           u &= g, && \text{on } \partial\Omega,
# \end{align*}
# where the domain is given using polar coordinates,
# \begin{gather*}
#   \Omega = \{ (r,\varphi)\colon r\in(0,1), \varphi\in(0,\Phi) \}~.
# \end{gather*}
# For the boundary condition $g$, we set it to the trace of the function $u$, given by
# \begin{gather*}
#   u(r,\varphi) = r^{\frac{\pi}{\Phi}} \sin\big(\frac{\pi}{\Phi} \varphi \big)
# \end{gather*}
#
# We first define the domain and set up the grid and space.

# <codecell>
try:
    get_ipython().magic(u'matplotlib inline # can also use notebook or nbagg')
except:
    pass
import math
import numpy
import matplotlib.pyplot as pyplot
from dune.fem.view import adaptiveLeafGridView
from dune.fem.plotting import plotPointData as plot
import dune.grid as grid
import dune.fem as fem
from dune.fem.view import adaptiveLeafGridView as gridView
from dune.fem.space import lagrange as solutionSpace
from dune.alugrid import aluConformGrid as hierachicalGrid


# set the angle for the corner (0<angle<=360)
cornerAngle = 320.

# use a second order space
order = 2

# <markdowncell>
# define the grid for this domain (vertices are the origin and 4
# equally spaced points on the unit sphere starting with (1,0) and
# ending at (cos(cornerAngle), sin(cornerAngle))

# <codecell>
vertices = numpy.zeros((8, 2))
vertices[0] = [0, 0]
for i in range(0, 7):
    vertices[i+1] = [math.cos(cornerAngle/6*math.pi/180*i),
                     math.sin(cornerAngle/6*math.pi/180*i)]
triangles = numpy.array([[2,1,0], [0,3,2], [4,3,0],
                         [0,5,4], [6,5,0], [0,7,6]])
domain = {"vertices": vertices, "simplices": triangles}
grid = gridView( hierachicalGrid(domain) )
grid.hierarchicalGrid.globalRefine(2)
space  = solutionSpace(grid, dimrange=1, order=order)

# <markdowncell>
# Next we define the model together with the exact solution.

# <codecell>
from ufl import *
from dune.ufl import DirichletBC
from dune.fem.scheme import galerkin as solutionScheme
u = TrialFunction(space)
v = TestFunction(space)
x = SpatialCoordinate(space.cell())

# exact solution for this angle
Phi = cornerAngle / 180 * pi
phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*pi, 0)
exact = as_vector([inner(x, x)**(pi/2/Phi) * sin(pi/Phi * phi)])
a = inner(grad(u), grad(v)) * dx

# set up the scheme
laplace = solutionScheme([a==0, DirichletBC(space, exact, 1)])
uh = space.interpolate(lambda x: [0], name="solution")


# <markdowncell>
# Theory tells us that
# \begin{align*}
#   \int_\Omega |\nabla(u-u_h)|^2 \leq \sum_K \eta_K,
# \end{align*}
# where on each element $K$ of the grid the local estimator is given by
# \begin{align*}
#   \eta_K = h_K^2 \int_K |\triangle u_h|^2 +
#     \frac{1}{2}\sum_{S\subset \partial K} h_S \int_S [\nabla u_h]^2.
# \end{align*}
# Here $[\cdot]$ is the jump in normal direction over the edges of the grid.
#
# We compute the elementwise indicator by defining a bilinear form
# \begin{align*}
#   \eta(u,v) = \int_\Omega h^2 |\triangle u_h|^2 v +
#     \int_{I_h} h_S [\nabla u_h]^2 \{v\},
# \end{align*}
# where $\{\cdot\}$ is the average over the cell edges and $[\cdot]$
# the jump. With $h$ and $h_S$ we denote local grid spacings and with
# $I_h$ the set of all facets in the grid.
# This bilinear form can be easily written in UFL and by using it to
# define a discrete operator $L$ from the second order Lagrange space into a space containing piecewise constant functions
# we have $L[u_h]|_{K} = \eta_K$.

# <codecell>
# energy error
h1error = inner(grad(uh - exact), grad(uh - exact))

# residual estimator
from dune.fem.space import finiteVolume as estimatorSpace
from dune.fem.operator import galerkin as estimatorOp

fvspace = estimatorSpace(grid, dimrange=1)
estimate = fvspace.interpolate([0], name="estimate")

hT = MaxCellEdgeLength(space.cell())
he = MaxFacetEdgeLength(space.cell())('+')
n = FacetNormal(space.cell())
estimator_ufl = hT**2 * (div(grad(u[0])))**2 * v[0] * dx +\
        he * inner(jump(grad(u[0])), n('+'))**2 * avg(v[0]) * dS
estimator = estimatorOp(estimator_ufl, space, fvspace)
# marking strategy (equidistribution)
tolerance = 0.1


# <markdowncell>
# Let us solve over a loop (solve,estimate mark) and plot the solutions side by side.

# <codecell>
fig = pyplot.figure(figsize=(10,10))
count = 0
while count < 20:
    laplace.solve(target=uh)
    if count%3 == 0:
        pyplot.show()
        pyplot.close('all')
        fig = pyplot.figure(figsize=(10,10))
    plot(uh, figure=(fig, 131+count%3), colorbar=False)
    # compute the actual error and the estimator
    error = math.sqrt(fem.function.integrate(grid, h1error, 5)[0])
    estimator(uh, estimate)
    eta = sum(estimate.dofVector)
    print(count, ": size=", grid.size(0), "estimate=", eta,
          "error=", error)
    if eta < tolerance:
        break
    if tolerance == 0.:
        grid.hierarchicalGrid.globalRefine(2)
        uh.interpolate([0])  # initial guess needed
    else:
        marked = fem.mark(estimate,tolerance/grid.size(0))
        fem.adapt([uh])
        fem.loadBalance([uh])
    laplace.solve( target=uh )
    count += 1
pyplot.show()
pyplot.close('all')

# <markdowncell>
# Let's have a look at the center of the domain:

# <codecell>
fig = pyplot.figure(figsize=(15,15))
plot(uh, figure=(fig, 131), xlim=(-0.5, 0.5),
     ylim=(-0.5, 0.5), colorbar={"shrink": 0.25})
plot(uh, figure=(fig, 132), xlim=(-0.25, 0.25),
     ylim=(-0.25, 0.25),colorbar={"shrink": 0.25})
plot(uh, figure=(fig, 133), xlim=(-0.125, 0.125),
     ylim=(-0.125, 0.125),colorbar={"shrink": 0.25})
pyplot.show()
pyplot.close('all')


# <markdowncell>
# Finally, let us have a look at the grid levels.

# <codecell>
from dune.fem.function import levelFunction
plot(levelFunction(grid), xlim=(-0.2,1), ylim=(-0.2,1))
