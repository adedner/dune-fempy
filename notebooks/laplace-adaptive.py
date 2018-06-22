
# coding: utf-8

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

# In[1]:


try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass
import math
import numpy
import matplotlib.pyplot as pyplot
import dune.create as create
from dune.grid import gridFunction
from dune.fem.view import adaptiveLeafGridView
from dune.fem.plotting import plotPointData as plot
import dune.grid as grid
import dune.fem as fem

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
view = create.view("adaptive", "ALUConform", domain)
view.hierarchicalGrid.globalRefine(2)
spc  = create.space( "lagrange", view, dimrange=1, order=order )

# Next define the model together with the exact solution.

# In[2]:


from ufl import *
from dune.ufl import DirichletBC
u = TrialFunction(spc)
v = TestFunction(spc)
x = SpatialCoordinate(spc.cell())

# exact solution for this angle
Phi = cornerAngle / 180 * math.pi
phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*math.pi, 0)
exact = as_vector([inner(x,x)**(math.pi/2/Phi) * sin(math.pi/Phi * phi)])
a = inner(grad(u), grad(v)) * dx

# set up the scheme
laplace = create.scheme("h1", spc, [a==0, DirichletBC(spc,exact,1)])
uh = spc.interpolate(lambda x: [0], name="solution")


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

# In[3]:


# energy error
h1error = inner(grad(uh - exact), grad(uh - exact))

# residual estimator
fvspc = create.space("finitevolume", view, dimrange=1, storage="istl")
estimate = fvspc.interpolate([0], name="estimate")

hT = MaxCellEdgeLength(spc.cell())
he = MaxFacetEdgeLength(spc.cell())('+')
n = FacetNormal(spc.cell())
estimator_ufl = hT**2 * (div(grad(u[0])))**2 * v[0] * dx                 + he * inner(jump(grad(u[0])), n('+'))**2 * avg(v[0]) * dS
estimator = create.operator("galerkin", estimator_ufl, spc, fvspc)

# marking strategy (equidistribution)
tolerance = 0.1
def mark(element):
    estLocal = estimate(element, element.geometry.referenceElement.center)
    return grid.Marker.refine if estLocal[0] > tolerance / gridSize else grid.Marker.keep


# In[4]:


# adaptive loop (solve, mark, estimate)
# fig = pyplot.figure(figsize=(10,10))
count = 0
while count < 20:
    laplace.solve(target=uh)
    if False: # count%3 == 0:
        pyplot.show()
        pyplot.close('all')
        fig = pyplot.figure(figsize=(10,10))
    # plot(uh,figure=(fig,131+count%3), colorbar=False)
    # compute the actual error and the estimator
    error = math.sqrt(fem.function.integrate(view, h1error, 5)[0])
    estimator(uh, estimate)
    eta = sum(estimate.dofVector)
    gridSize = view.size(0)
    print(count, ": size=", gridSize, "estimate=", eta, "error=", error)
    if eta < tolerance:
        break
    if tolerance == 0.:
        view.hierarchicalGrid.globalRefine(2)
        uh.interpolate([0])  # initial guess needed
    else:
        marked = view.hierarchicalGrid.mark(mark)
        fem.adapt(view.hierarchicalGrid, [uh])
        fem.loadBalance(view.hierarchicalGrid, [uh])
    count += 1
# pyplot.show()
# pyplot.close('all')


# Let's have a look at the center of the domain:

# In[5]:


fig = pyplot.figure(figsize=(15,15))
plot(uh, figure=(fig,131+0), xlim=(-0.5,0.5), ylim=(-0.5,0.5),colorbar={"shrink":0.3})
plot(uh, figure=(fig,131+1), xlim=(-0.25,0.25), ylim=(-0.25,0.25),colorbar={"shrink":0.3})
plot(uh, figure=(fig,131+2), xlim=(-0.125,0.125), ylim=(-0.125,0.125),colorbar={"shrink":0.3})
pyplot.show()
pyplot.close('all')


# Finally, let us have a look at the grid levels:

# In[ ]:


from dune.fem.function import levelFunction
# plot(levelFunction(view), xlim=(-0.2,1), ylim=(-0.2,1))


# ## HP-Adaptivity
#
# At the moment `dune-fempy` only supports p-adaptivity for discontinuous spaces. So first we need to formulate the problem in a dg context and adapt the residual estimator accordingly:

# In[ ]:

maxOrder = 4
# view = create.view("adaptive", "ALUConform", domain, dimgrid=2)
# dgSpc  = create.space( "dgonbhp", view, dimrange=1, order=maxOrder )
# make a remark on vertex order
# make a remark on non affine cubes
vertices = numpy.array([(0,0), (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1), (1,-1), (1,0)])
cubes = numpy.array([[0,1,3,2], [0,3,5,4], [0,5,7,6], [0,9,7,8]])
domain = {"vertices": vertices, "cubes": cubes}
view = create.view("adaptive", "ALUCube", domain)
plot(view)
dgSpc  = create.space( "dglegendrehp", view, dimrange=1, order=maxOrder )
view.hierarchicalGrid.globalRefine(2)

# exact solution for this angle
cornerAngle = 360.
u = TrialFunction(dgSpc)
v = TestFunction(dgSpc)
x = SpatialCoordinate(dgSpc.cell())
Phi = cornerAngle / 180 * math.pi
phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*math.pi, 0)
exact = as_vector([inner(x,x)**(math.pi/2/Phi) * sin(math.pi/Phi * phi)])
n, h = FacetNormal(dgSpc.cell()), MinFacetEdgeLength(dgSpc.cell())
mu   = 100. / avg(h)
a  = inner(grad(u), grad(v)) * dx
a -= inner(outer(jump(u), n('+')), avg(grad(v))) * dS
a -= inner(avg(grad(u)), outer(jump(v), n('+'))) * dS
a += mu * inner(jump(u), jump(v)) * dS
# a -= inner(outer(u, n), grad(v)) * ds
# a -= inner(grad(u), outer(v, n)) * ds
a += mu * inner(u-exact, v) * ds
a += div(grad(exact[0])) * v[0] * dx

# laplace = create.scheme("dg", dgSpc, [a==0, DirichletBC(spc,exact,1)], penalty=20)
laplace = create.scheme("galerkin", dgSpc, a==0)
uh = dgSpc.interpolate([0],name="solution")
uh_pm1 = dgSpc.interpolate([0],name="solution_pm1")

# energy error
h1error = inner(grad(uh - exact), grad(uh - exact))

# residual estimator
fvspc = create.space("finitevolume", view, dimrange=1, storage="istl")
estimate = fvspc.interpolate([0], name="estimate")
estimate_pm1 = fvspc.interpolate([0], name="estimate_pm1")

hT = MaxCellEdgeLength(dgSpc.cell())
he = MaxFacetEdgeLength(dgSpc.cell())('+')
estimator_ufl = hT**2 * (div(grad(u[0])))**2 * v[0] * dx +\
                he * inner(jump(grad(u[0])), n('+'))**2 * avg(v[0]) * dS +\
                1/he * jump(u[0])**2 * dS +\
                1/he * (u[0]-exact[0])**2 * ds
estimator = create.operator("galerkin", estimator_ufl, dgSpc, fvspc)

# Note: there are not enough test here
# 1: one can use spaceAdapt with a space that is not padaptive and does not fit the df vector (e.g. spc)
# 2: compiler error if the OrderReduction is not used with a p-adaptive space
from dune.femdg import createOrderRedcution
orderreduce = createOrderRedcution( dgSpc )
fem.spaceAdapt(dgSpc, lambda e: 2, [uh])

# In[ ]:
pTol = 0.01
orders = [0,]*(maxOrder+1)
def markp(element):
    return 2
    polorder = dgSpc.localOrder(element)
    center = element.geometry.referenceElement.center
    r      = polorder*estimate(element,center)[0]
    r_p1   = (polorder-1)*estimate_pm1(element,center)[0] * element.geometry.volume
    if r_p1 < 0.1*pTol*r:
        ret = polorder-1 if polorder > 2 else polorder
    elif r_p1+1e-15 >= pTol*r:
        ret = polorder+1 if polorder < maxOrder else polorder
    else:
        ret = polorder
    orders[ret] += 1
    return ret
pMarker = gridFunction(view)(lambda e,x: markp(e))

# adaptive loop (solve, mark, estimate)
fig = pyplot.figure(figsize=(10,10))
count = 0
tolerance = 0.05
while count < 20:
    print("solve",flush=True)
    laplace.solve(target=uh)
    plot(uh,figure=(fig,141), colorbar=False)
    # compute the actual error and the estimator
    error = math.sqrt(fem.function.integrate(view, h1error, 5)[0])
    estimator(uh, estimate)
    eta = sum(estimate.dofVector)
    gridSize = view.size(0)
    print(count, ": size=", gridSize, "estimate=", eta, "error=", error)
    if eta < tolerance:
        break
    if tolerance == 0.:
        view.hierarchicalGrid.globalRefine(2)
        uh.interpolate([0])  # initial guess needed
    else:
        print("h-adapt",flush=True)
        marked = view.hierarchicalGrid.mark(mark)
        fem.adapt(view.hierarchicalGrid, [uh])
        fem.loadBalance(view.hierarchicalGrid, [uh])
        print("p-adapt",flush=True)
        estimator(uh, estimate)
        orderreduce(uh,uh_pm1)
        estimator(uh_pm1, estimate_pm1)
        logScale = lambda u: gridFunction(view)(lambda e,x: math.log10(u(e,x)))
        plot(logScale(estimate),figure=(fig,142), colorbar=False)
        plot(logScale(estimate_pm1),figure=(fig,143), colorbar=False)
        plot(pMarker,figure=(fig,144), colorbar=True)
        for i in range(len(orders)): orders[i] = 0
        fem.spaceAdapt(dgSpc, markp, [uh])
        print(orders)

    count += 1
    pyplot.show()
    pyplot.close('all')
    fig = pyplot.figure(figsize=(10,10))
pyplot.show()
pyplot.close('all')
