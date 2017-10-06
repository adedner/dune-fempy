# coding: utf-8

# # Heat Equation - adding coefficient and constants to the model [(Notebook)][1]
#
# [1]: _downloads/laplace-coefficients.ipynb
#

# In[1]:

from __future__ import print_function
try:
    get_ipython().magic(u'matplotlib inline # can also use notebook or nbagg')
except:
    pass

import math
from ufl import *

from dune.grid import cartesianDomain

import dune.create as create
import dune.fem
from dune.fem.plotting import plotPointData as plot

dune.fem.parameter.append({"fem.verboserank": 0,
                           "fem.preconditioning.method": "ilu-0",
                           "fem.preconditioning.iterations": 1,
                           "fem.preconditioning.relaxation": 1.2})


# In[2]:

# Crank Nicholson
theta = 0.5

# set up a 2d simplex grid over the interval [0,1]^2 with h = 1/16
grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
# set up a lagrange scalar space with polynomial order 2 over that grid
spc = create.space("Lagrange", grid, dimrange=1, order=2, storage="fem")

# set up initial conditions
solution = spc.interpolate(lambda x: [math.atan((10.0 * x[0] * (1-x[0]) * x[1] * (1-x[1]))**2)], name="u")
plot(solution)


# In[3]:

# get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
old_solution = solution.copy();

# now define the actual pde to solve:
#            u - u_n deltaT laplace( theta u + (1-theta) u_n ) = 0
u = spc.uflTrialFunction()
v = spc.uflTestFunction()
tau = spc.uflNamedConstant(name="tau")
a = (inner(u - old_solution, v) +    tau * inner(grad(theta*u + (1-theta)*old_solution), grad(v)) ) * dx

# now generate the model code and compile
model = create.model("elliptic", grid, a == 0)

# setup structure for olver parameters
solverParameter={"fem.solver.newton.linabstol": 1e-13,
                 "fem.solver.newton.linreduction": 1e-13,
                 "fem.solver.newton.tolerance": 1e-12,
                 "fem.solver.newton.verbose": "true",
                 "fem.solver.newton.linear.verbose": "false"}
# create the solver using a standard fem scheme
scheme = create.scheme("h1", spc, model, parameters=solverParameter)


# In[4]:

endTime = 0.4
deltaT = 0.01
model.tau = deltaT

# now loop through time and output the solution after each time step
steps = int(endTime / deltaT)
for n in range(1,steps+1):
    old_solution.assign(solution)
    scheme.solve(target=solution)
    if n % 4 == 3:
        plot(solution)
    # grid.writeVTK("heat", pointdata=[solution], number=n)
