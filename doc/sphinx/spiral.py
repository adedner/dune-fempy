# coding: utf-8

# #  Spiral Wave

# In[1]:

from __future__ import print_function
import math
from functools import reduce

import ufl

import dune.ufl
import dune.grid
import dune.models.elliptic
import dune.fem
import dune.create as create

# http://www.scholarpedia.org/article/Barkley_model
dimRange   = 2
dt         = 0.25
if 1:
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

def initial(x):
    return [ 1   if x[1]>1.25 else 0,             0.5 if x[0]<1.25 else 0. ]


# In[2]:

# set up reference domain
domain = dune.grid.cartesianDomain([0,0],[3.5,3.5],[40,40])
domain = dune.grid.cartesianDomain([0,0],[2.5,2.5],[30,30])
grid = create.grid("Yasp", domain, dimgrid=2)
spc  = create.space( "Lagrange", grid, dimrange=dimRange, order=1 )

solution    = spc.interpolate( initial, name="solution" )
solution_n  = spc.interpolate( initial, name="solution_n" )
forcing     = spc.interpolate( [0,0,0], name="forcing" )


# In[3]:

uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), dimRange)
u = ufl.TrialFunction(uflSpace)
v = ufl.TestFunction(uflSpace)
un = dune.ufl.NamedCoefficient(uflSpace, "un")

# right hand sie (time derivative part + explicit forcing in v)
a_ex = (ufl.inner(un, v) + dt * ufl.inner(spiral_h(un[0], un[1]), v[1])) * ufl.dx
# left hand side (heat equation in first variable + backward Euler in time)
a_im = (dt * spiral_D * ufl.inner(ufl.grad(u[0]), ufl.grad(v[0])) + ufl.inner(u,v)) * ufl.dx

modelCode = dune.models.elliptic.compileUFL(a_im == a_ex)

# extra model source code
# -----------------------
sourceCode = """      double uth = (@gf:un[ 1 ] + @const:b) / @const:a;
      if( @gf:un[ 0 ] <= uth )
        result[ 0 ] -= @const:dt/@const:eps * u[ 0 ] * (1.0 - @gf:un[ 0 ]) * (@gf:un[ 0 ] - uth);
      else
        result[ 0 ] -= @const:dt/@const:eps * @gf:un[ 0 ] * (1.0 - u[ 0 ]) * (@gf:un[ 0 ] - uth);
"""
linSourceCode = """      double uth = (@gf:un[ 1 ] + @const:b) / @const:a;
      if( @gf:un[ 0 ] <= uth )
        result[ 0 ] -= @const:dt/@const:eps * u[ 0 ] * (1.0 - @gf:un[ 0 ]) * (@gf:un[ 0 ] - uth);
      else
        result[ 0 ] -= @const:dt/@const:eps * @gf:un[ 0 ] * (-u[ 0 ]) * (@gf:un[ 0 ] - uth);
"""
modelCode.appendCode('source', sourceCode, coefficients={"un": solution_n} )
modelCode.appendCode('linSource', linSourceCode )


# In[4]:

model = create.model("elliptic", grid, modelCode, coefficients={"un": solution_n} )
model.setConstant("a", [spiral_a])
model.setConstant("b", [spiral_b])
model.setConstant("eps", [spiral_eps])
model.setConstant("dt", [dt])
solverParameters = {
        "fem.solver.newton.tolerance": 1e-3,
        "fem.solver.newton.linabstol": 1e-5,
        "fem.solver.newton.linreduction": 1e-5,
        "fem.solver.newton.verbose": 0,
        "fem.solver.newton.linear.verbose": 0
    }
scheme = create.scheme("h1", spc, model, ("pardg","cg"),parameters=solverParameters)


# In[5]:

import matplotlib.pyplot as plt
from numpy import linspace
from matplotlib import animation, rc
rc('animation', html='html5')
fig, ax = plt.subplots()
ax.set_xlim(( 0, 2.5))
ax.set_ylim(( 0, 2.5))
triangulation = grid.triangulation(1)
levels = linspace(-0.1, 1.1, 256)
ax.set_aspect('equal')
t        = 0.
stepsize = 0.5
nextstep = 0.
iterations = 0

def animate(count):
    global t, stepsize, nextstep, iterations, dt, solution
    # print('frame',count,t)
    while t < nextstep:
        solution_n.assign(solution)
        _,info = scheme.solve(target=solution)
        # print("Computing solution a t = " + str(t + dt), "iterations: " + info["linear_iterations"] )
        iterations += int( info["linear_iterations"] )
        t     += dt
    data = solution.pointData(1)
    plt.tricontourf(triangulation, data[:,0], cmap=plt.cm.rainbow, levels=levels)
    # grid.writeVTK("spiral", pointdata=[solution], number=count)
    nextstep += stepsize


# In[6]:

animation.FuncAnimation(fig, animate, frames=25, interval=100, blit=False)
