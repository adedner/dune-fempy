# coding: utf-8

# # Adaptive Finite Element
# We study the classic _reentrand corner_ problem:
# \begin{align*}
#   -\Delta u &= f && \text{in } \Omega \\
#           u &= g && \text{on } \partial\Omega
# \end{align*}
# where the domain is given using polar coordinates
# \begin{gather}
#   \Omega = \{ (r,\varphi)\colon r\in(0,1) \varphi\in(0,\Phi) \}~.
# \end{gather}
# For $g$ we take the trace of the function u, given by
# \begin{gather}
#   u(r,\varphi) = r^{\frac{\pi}{\Phi}} \sin\big(\frac{\pi}{\Phi} \varphi \big)
# \end{gather}
#
# We first define the domain and set up the grid and space

# In[1]:

import math
import dune.create as create
from dune.fem.view import adaptiveLeafGridView
import dune.grid as grid
import dune.fem as fem

# set the angle for the corner (0<angle<=360)
cornerAngle = 320.

# use a second order space
order = 2

# define the grid for this domain (vertices are the origin and 4 equally spaces on the
# unit sphere starting with (1,0) and ending at # (cos(cornerAngle),sin(cornerAngle))
vertices = [(0,0)]
dgf = "VERTEX\n 0 0 \n"
for i in range(0,5):
    dgf += str(math.cos(cornerAngle/4*math.pi/180*i)) + " " + str(math.sin(cornerAngle/4*math.pi/180*i)) + "\n"
dgf += "#\n"
dgf += """
SIMPLEX
  0 1 2
  0 2 3
  0 3 4
  0 4 5
#
BOUNDARYDOMAIN
  default 1
#
PROJECTION
  function p(x) = x / |x|
  segment  1 2  p
  segment  2 3  p
  segment  3 4  p
  segment  4 5  p
#
"""
grid = create.view("adaptive", "ALUConform", grid.string2dgf(dgf), dimgrid=2)
grid.hierarchicalGrid.globalRefine(2)
spc  = create.space( "Lagrange", grid, dimrange=1, order=order )


# Next define the model together with the exact solution

# In[3]:

from ufl import *
from dune.ufl import Space
uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

# exact solution for this angle
def exact(x):
    Phi = cornerAngle / 180 * math.pi
    r2 = x.two_norm2
    phi = math.atan2(x[1], x[0])
    if x[1] < 0:
        phi += 2*math.pi
    return [r2**(math.pi/2/Phi) * sin(math.pi/Phi * phi)]
def exactJac(x):
    r2 = x.two_norm2
    phi = math.atan2(x[1], x[0])
    if x[1] < 0:
        phi += 2*math.pi
    r2dx=2.*x[0]
    r2dy=2.*x[1]
    if r2 == 0:
        phidx = 0
        phidy = 0
        r2pow = 0
    else:
        phidx=-x[1]/r2
        phidy=x[0]/r2
        r2pow = r2**(90./cornerAngle-1)
    dx = r2pow * ( 90./cornerAngle*r2dx * sin(180./cornerAngle * phi)
                  + r2 * 180./cornerAngle*cos( 180./cornerAngle * phi) * phidx )
    dy = r2pow * ( 90./cornerAngle*r2dy * sin(180./cornerAngle * phi)
                  + r2 * 180./cornerAngle*cos( 180./cornerAngle * phi) * phidy )
    return [dx,dy]

exact_gf = create.function("global", grid, "exact", order+1, exact)
bnd_u = Coefficient(uflSpace)
a = inner(grad(u), grad(v)) * dx
model = create.model("elliptic", grid, a == 0, dirichlet={1: bnd_u}, tempVars=False, coefficients={bnd_u: exact_gf})


# In[5]:

try:
    import matplotlib
    from matplotlib import pyplot
    from numpy import amin, amax, linspace
    from IPython.core.display import display

    def plot(grid, solution, xlim = None, ylim = None):
        triangulation = grid.triangulation(4)
        data = solution.pointData(4)

        levels = linspace(amin(data[:,0]), amax(data[:,0]), 256)

        fig = pyplot.figure()
        fig.gca().set_aspect('equal')
        if xlim:
            fig.gca().set_xlim(xlim)
        if ylim:
            fig.gca().set_ylim(ylim)
        pyplot.triplot(grid.triangulation(), antialiased=True, linewidth=0.2, color='black')
        pyplot.tricontourf(triangulation, data[:,0], cmap=pyplot.cm.rainbow, levels=levels)
        display(pyplot.gcf())
except ImportError as e:
    print(e)
    def plot(grid, solution):
        pass


# In[6]:

pyplot.close("all")
# set up the scheme
laplace = create.scheme("h1", spc, model)
uh = spc.interpolate(lambda x: [0], name="solution")
laplace.solve(target=uh)

# function used for computing approximation error
def h1error(en,x):
    y = en.geometry.position(x)
    val = uh.localFunction(en).evaluate(x) - exact(y)
    jac = uh.localFunction(en).jacobian(x)[0] - exactJac(y)
    return [ sqrt( val[0]*val[0] + jac*jac) ];
h1error_gf = create.function( "local", grid, "error", order+1, h1error )

# adaptive loop (mark, estimate, solve)
count = 0
tol = 0.05 # use 0 for global refinement
while count < 8:
    error = math.sqrt(h1error_gf.integrate()[0])
    [estimate, marked] = laplace.mark(uh, tol)
    plot(grid, uh)
    print(count, ": size=",grid.size(0), "estimate=",estimate,"error=",error)
    if marked == False or estimate < tol:
        break
    if tol == 0.:
        grid.hierarchicalGrid.globalRefine(2)
        uh.interpolate([0])  # initial guess needed
    else:
        fem.adapt(grid.hierarchicalGrid, [uh])
        fem.loadBalance(grid.hierarchicalGrid, [uh])
    laplace.solve( target=uh )
    count += 1


# Let's have a look at the center of the domain:

# In[7]:

plot(grid,uh, (-0.5,0.5),(-0.5,0.5))
plot(grid,uh, (-0.25,0.25),(-0.25,0.25))
plot(grid,uh, (-0.125,0.125),(-0.125,0.125))
plot(grid,uh, (-0.0625,0.0625),(-0.0625,0.0625))


# In[ ]:
