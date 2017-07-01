# coding: utf-8

# # DG Schemes [(Notebook)][1]
#
# [1]: _downloads/laplace-dg.ipynb
# show the
# - dgscheme
# - galerkin schemes

# Let us consider a simple Laplace problem with Dirichlet boundary conditions:
# \begin{equation*}
#   \begin{aligned}
#     -\Delta u &= \sin(\pi x_1) \sin(\pi x_2) && \text{in $\Omega$}, \\
#             u &= 0 && \text{on $\partial\Omega$}.
#   \end{aligned}
# \end{equation*}
#
# First, we need to set up a computational grid and a discontinuous ansatz space on it. Here, we use the orthonormal discontinuous space:

# In[1]:

try:
    get_ipython().magic(u'matplotlib inline # can also use notebook or nbagg')
except:
    pass
from dune.grid import cartesianDomain
from dune.fem import parameter
from dune.fem.plotting import plotPointData as plot

import dune.create as create

parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu-0", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
spc = create.space("DGONB", grid, dimrange=1, order=2, storage="istl")


# The classical IPDG method for this problem reads
# \begin{equation*}
#   \int_\Omega \nabla u\,\nabla \varphi\,dx
#     - \int_\Gamma ([[u]] \otimes \vec{n} : \{\{\nabla \varphi\}\} + \{\{\nabla u\}\} : [[\varphi]] \otimes \vec{n}\,dx
#     + \int_\Gamma \frac{\mu}{h} [[u]] [[\varphi]]
#     = 0.
# \end{equation*}
#
# The following code implements this equation in UFL notation:

# In[2]:

import math
from ufl import *

from dune.ufl import Space

uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
n, h = FacetNormal(uflSpace.cell()), MinFacetEdgeLength(uflSpace.cell())
mu = 7.5 / avg(h)

a = inner(grad(u), grad(v)) * dx
a -= (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
a += mu * inner(jump(u), jump(v)) * dS
a -= (inner(outer(u, n), grad(v)) + inner(grad(u), outer(v, n))) * ds
a += mu * inner(u, v) * ds

b = sin(pi*x[0])*sin(pi*x[1])*v[0]*dx


# Next, we compile this into the *integrands*, plug them into the *galerkin* scheme and solve the problem:

# In[3]:

model = create.model("integrands", grid, a == b)

newtonParameter = {"linabstol": 1e-13, "linreduction": 1e-13, "tolerance": 1e-12, "verbose": "true", "linear.verbose": "false"}
scheme = create.scheme("galerkin", spc, model, parameters={"fem.solver.newton." + k: v for k, v in newtonParameter.items()})

uh, _ = scheme.solve()


# The result looks as follows:

# In[4]:

plot(uh)
