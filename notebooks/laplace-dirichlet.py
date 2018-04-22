
# coding: utf-8

# # Dirichlet Boundary Conditions [(Notebook)][1]
#
# [1]: _downloads/laplace-dirichlet.ipynb
#
# __Still need to be able to set individual boundary ids....__
#
# Consider a simple elliptic problem with Dirichlet and Neuman boundary conditions
# \begin{align*}
#   -\Delta u &= f && \text{in } \Omega \\
#           u &= g && \text{on } \Gamma_D \\
#            -\nabla u\cdot n + \alpha u &= r && \text{on } \Gamma_R
# \end{align*}
# Note that $\alpha$ can be zero so that this includes non zero Neuman boundary conditions.
#
# `dune-fempy` uses integer boundary identifies to distinguish between different parts of the boundary. These are set during grid construction on the coarse grid boundaries:

# In[ ]:


from __future__ import print_function, division

try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass

import math
import numpy

from ufl import *

from dune.grid import cartesianDomain
from dune.fem import parameter
from dune.fem.plotting import plotPointData as plot
from dune.ufl import Space, DirichletBC

import dune.create as create

parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu-0", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})


# In[ ]:


vertices = numpy.array([(0,0), (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1)])
triangles = numpy.array([(0,1,2), (0,2,3), (0,3,4), (0,4,5), (0,5,6), (0,6,7)])

grid = create.grid("ALUConform", {"vertices": vertices, "simplices": triangles}, dimgrid=2)
grid.hierarchicalGrid.globalRefine(4)


# In[ ]:


spc = create.space("lagrange", grid, dimrange=1, order=1, storage="istl")

uflSpace = Space(spc)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*math.pi, 0)
exact = as_vector([inner(x,x)**(0.5*180/270) * sin((180/270) * phi)])
a = inner(grad(u), grad(v))*dx


# In[ ]:


newtonParameter = {"linabstol": 1e-13, "linreduction": 1e-13, "tolerance": 1e-12, "verbose": "true", "linear.verbose": "false"}
scheme = create.scheme("h1", spc, [a==0,DirichletBC(uflSpace,exact,1)],
                parameters={"fem.solver.newton." + k: v for k, v in newtonParameter.items()})

solution, _ = scheme.solve()
plot(solution)
