
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

# In[ ]:


try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass
from dune.grid import cartesianDomain
from dune.fem import parameter
from dune.fem.function import integrate
from dune.fem.plotting import plotPointData as plot

import dune.create as create

parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})
order = 2

grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
spc = create.space("dgonb", grid, dimrange=1, order=order, storage="istl")


# The classical IPDG method for this problem reads
# \begin{equation*}
#   \int_\Omega \nabla u\,\nabla \varphi\,dx
#     - \int_\Gamma ([[u]] \otimes \vec{n} : \{\{\nabla \varphi\}\} + \{\{\nabla u\}\} : [[\varphi]] \otimes \vec{n}\,dx
#     + \int_\Gamma \frac{\mu}{h} [[u]] [[\varphi]]
#     = 0.
# \end{equation*}
#
# The following code implements this equation in UFL notation:

# In[ ]:


import math
from ufl import *

u = TrialFunction(spc)
v = TestFunction(spc)
x = SpatialCoordinate(spc)
n, h = FacetNormal(spc), MinFacetEdgeLength(spc)
mu = 7.5 / avg(h)

a = inner(grad(u), grad(v)) * dx
a -= (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
a += mu * inner(jump(u), jump(v)) * dS
a -= (inner(outer(u, n), grad(v)) + inner(grad(u), outer(v, n))) * ds
a += mu * inner(u, v) * ds

b = 2*pi*pi*sin(pi*x[0])*sin(pi*x[1])*v[0]*dx
exact = as_vector([ sin(pi*x[0])*sin(pi*x[1]) ])

# Next, we compile this into the *integrands*, plug them into the *galerkin* scheme and solve the problem:

# In[ ]:


newtonParameter = {"tolerance": 1e-5, "verbose": "false",
                   "linear.absolutetol": 1e-8, "linear.reductiontol": 1e-8,
                   "linear.preconditioning.method": "ilu",
                   "linear.preconditioning.iterations": 1, "linear.preconditioning.relaxation": 1.2,
                   "linear.verbose": "false"}
scheme = create.scheme("galerkin", a == b, spc, parameters={"newton." + k: v for k, v in newtonParameter.items()})

uh = spc.interpolate([0],name="dg")
scheme.solve(target=uh)


# The result looks as follows:

# In[ ]:


plot(uh)
error = uh - exact
print("DG: L^2 and H^1 error:",
  [ sqrt(e) for e in integrate(grid,[error**2,inner(grad(error),grad(error))], order=5) ] )

# We can simply produce a continuous projection of the discontinuous function
lagSpc = create.space("lagrange", grid, dimrange=1, order=order+1, storage="istl")
contUh = lagSpc.project(uh, name="Vtx")
plot(contUh)
error = contUh - exact
print("Lag: L^2 and H^1 error:",
  [ sqrt(e) for e in integrate(grid,[error**2,inner(grad(error),grad(error))], order=5) ] )

# using averaged gradients
from dune.ufl import DirichletBC
gradSpc = create.space("lagrange", grid, dimrange=grid.dimension, order=order, storage="istl")
contGrad = gradSpc.project(grad(uh[0]),name="ZZ")
a = inner(grad(u[0])-contGrad, grad(v[0])) * dx
scheme = create.scheme("galerkin", [a==0, DirichletBC(lagSpc,uh,1)],
    lagSpc, parameters={"newton." + k: v for k, v in newtonParameter.items()})

zzUh = lagSpc.interpolate([0,]*grid.dimension,name="ZZ")
scheme.solve(target=zzUh)
plot(zzUh)
error = zzUh - exact
print("ZZ: L^2 and H^1 error:",
  [ sqrt(e) for e in integrate(grid,[error**2,inner(grad(error),grad(error))], order=5) ] )

grid.writeVTK("reconstruction", pointdata=[uh,contUh,zzUh])
