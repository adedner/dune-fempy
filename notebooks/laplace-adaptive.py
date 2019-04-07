# coding: utf-8

# <markdowncell>
# # Adaptive Finite Element [(Notebook)][1]
#
# [1]: _downloads/laplace-adaptive.ipynb
# We study the classic _re-entrant corner_ problem:
# \begin{align*}
#   -\Delta u &= f && \text{in } \Omega \\
#           u &= g && \text{on } \partial\Omega
# \end{align*}
# where the domain is given using polar coordinates
# \begin{gather}
#   \Omega = \{ (r,\varphi)\colon r\in(0,1), \varphi\in(0,\Phi) \}~.
# \end{gather}
# For $g$ we take the trace of the function $u$, given by
# \begin{gather}
#   u(r,\varphi) = r^{\frac{\pi}{\Phi}} \sin\big(\frac{\pi}{\Phi} \varphi \big)
# \end{gather}
#
# We first define the domain and set up the grid and space
# <codecell>


try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass
import math
import numpy
import matplotlib.pyplot as pyplot
from dune.fem.view     import adaptiveLeafGridView as gridView
from dune.fem.space    import lagrange      as solutionSpace
from dune.fem.space    import finiteVolume  as estimatorSpace
from dune.fem.scheme   import galerkin      as solutionScheme
from dune.fem.operator import galerkin      as estimatorOperator
from dune.fem.plotting import plotPointData as plot
from dune.plotting     import block
import dune.grid as grid
import dune.fem  as fem

# set the angle for the corner (0<angle<=360)
cornerAngle = 320.

# use a second order space
order = 2

# define the grid for this domain (vertices are the origin and 4 equally spaced points on the
# unit sphere starting with (1,0) and ending at # (cos(cornerAngle), sin(cornerAngle))
vertices = numpy.zeros((8,2))
vertices[0] = [0,0]
for i in range(0,7):
    vertices[i+1] = [math.cos(cornerAngle/6*math.pi/180*i),
                     math.sin(cornerAngle/6*math.pi/180*i)]
triangles = numpy.array([[2,1,0], [0,3,2], [4,3,0], [0,5,4], [6,5,0], [0,7,6]])
domain = {"vertices": vertices, "simplices": triangles}
view = gridView("ALUConform", domain)
view.hierarchicalGrid.globalRefine(2)
space  = solutionSpace( view, dimRange=1, order=order, storage="fem" )

# <markdowncell>
# Next define the model together with the exact solution.
# <codecell>


from ufl import *
from dune.ufl import DirichletBC
u = TrialFunction(space)
v = TestFunction(space)
x = SpatialCoordinate(space.cell())

# exact solution for this angle
Phi = cornerAngle / 180 * math.pi
phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*math.pi, 0)
exact = as_vector([inner(x,x)**(math.pi/2/Phi) * sin(math.pi/Phi * phi)])
a = inner(grad(u), grad(v)) * dx

# set up the scheme
laplace = solutionScheme([a==0, DirichletBC(space,exact,1)], space)
uh = space.interpolate([0], name="solution")

# <markdowncell>
# Theory tells us that
# \begin{align*}
#   \int_\Omega \nabla(u-u_h) \leq \sum_K \eta_K
# \end{align*}
# where on each element $K$ of the grid the local estimator is given by
# \begin{align*}
#   \eta_K = h_K^2 \int_K |\triangle u_h|^2 +
#     \frac{1}{2}\sum_{S\subset \partial K} h_S \int_S [\nabla u_h]^2
# \end{align*}
# Here $[\cdot]$ is the jump in normal direction over the edges of the grid.
#
# We compute the elementwise indicator by defining a bilinear form
# \begin{align*}
#   \eta(u,v) = h_K^2 \int_K |\triangle u_h|^2 v +
#     \sum_{S\subset \partial K} h_S \int_S [\nabla u_h]^2 \{v\}
# \end{align*}
# where $\{\cdot\}$ is the average over the cell edges. This bilinear form can be easily written in UFL and by using it to define a discrete operator $L$ from the second order Lagrange space into a space containing piecewise constant functions
# we have $L[u_h]|_{K} = \eta_K$.
# <codecell>

# energy error
h1error = inner(grad(uh - exact), grad(uh - exact))

# residual estimator
fvspace = estimatorSpace(view, dimRange=1)
estimate = fvspace.interpolate([0], name="estimate")

hT = MaxCellEdgeLength(space.cell())
he = MaxFacetEdgeLength(space.cell())('+')
n = FacetNormal(space.cell())
estimator_ufl = hT**2 * (div(grad(u[0])))**2 * v[0] * dx + he * inner(jump(grad(u[0])), n('+'))**2 * avg(v[0]) * dS
estimator = estimatorOperator(estimator_ufl, space, fvspace)

# marking strategy (equidistribution)
tolerance = 0.1
gridSize = view.size(0)
def mark(element):
    estLocal = estimate(element, element.geometry.referenceElement.center)
    return grid.Marker.refine if estLocal[0] > tolerance / gridSize else grid.Marker.keep

# <markdowncell>
# Let's look at the result:
# <codecell>


# adaptive loop (solve, mark, estimate)
fig = pyplot.figure(figsize=(10,10))
count = 0
while count < 20:
    laplace.solve(target=uh)
    if count%3 == 0:
        pyplot.show(block=block)
        pyplot.close('all')
        fig = pyplot.figure(figsize=(10,10))
    plot(uh,figure=(fig,131+count%3), colorbar=False)
    # compute the actual error and the estimator
    error = math.sqrt(fem.function.integrate(view, h1error, 5))
    estimator(uh, estimate)
    eta = sum(estimate.dofVector)
    print(count, ": size=", gridSize, "estimate=", eta, "error=", error)
    if eta < tolerance:
        break
    if tolerance == 0.:
        view.hierarchicalGrid.globalRefine(2)
        uh.interpolate([0])  # initial guess needed
    else:
        marked = view.hierarchicalGrid.mark(mark)
        fem.adapt(uh)
        fem.loadBalance(uh)
    gridSize = view.size(0)
    count += 1

pyplot.show(block=block)
pyplot.close('all')

# <markdowncell>
# Let's have a look at the center of the domain:
# <codecell>


fig = pyplot.figure(figsize=(15,15))
plot(uh, figure=(fig,131+0), xlim=(-0.5,0.5), ylim=(-0.5,0.5),colorbar={"shrink":0.3})
plot(uh, figure=(fig,131+1), xlim=(-0.25,0.25), ylim=(-0.25,0.25),colorbar={"shrink":0.3})
plot(uh, figure=(fig,131+2), xlim=(-0.125,0.125), ylim=(-0.125,0.125),colorbar={"shrink":0.3})
pyplot.show(block=block)
pyplot.close('all')

# <markdowncell>
# Finally, let us have a look at the grid levels:
# <codecell>

from dune.fem.function import levelFunction
plot(levelFunction(view), xlim=(-0.2,1), ylim=(-0.2,1))
