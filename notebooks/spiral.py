
# coding: utf-8

# # Spiral Wave [(Notebook)][1]
#
# [1]: _downloads/spiral.ipynb

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

# Let's get started by importing some standard python packages, ufl, and some part of the dune-fempy package:

# In[ ]:


from __future__ import print_function
try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass
import math
from functools import reduce

import ufl

import dune.ufl
import dune.grid
import dune.models.integrands
import dune.fem
import dune.create as create


# In our attempt we will discretize the model as a 2x2 system. Here are some possible model parameters and initial conditions (we even have two sets of model parameters to choose from):

# In[ ]:


dimRange   = 1
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

initial_u = lambda x: [1   if x[1]>1.25 else 0]
initial_v = lambda x: [0.5 if x[0]<1.25 else 0]


# Now we set up the reference domain, the Lagrange finite element space (second order), and discrete functions for $(u^n,v^n($, $(u^{n+1},v^{n+1})$:

# In[ ]:


# domain = dune.grid.cartesianDomain([0,0],[3.5,3.5],[40,40])
domain = dune.grid.cartesianDomain([0,0],[2.5,2.5],[30,30])
grid = create.grid("ALUCube", domain, dimgrid=2)
spc  = create.space( "lagrange", grid, dimrange=dimRange, order=1 )

uh   = spc.interpolate( initial_u, name="u" )
uh_n = uh.copy()
vh   = spc.interpolate( initial_v, name="v" )
vh_n = vh.copy()


# We define the model in two steps:
# - first we define the standard parts, not involving $f_E,f_I$:
# - then we add the missing parts with the required _if_ statement directly using C++ code

# In[ ]:


uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), dimRange)
u   = ufl.TrialFunction(uflSpace)
phi = ufl.TestFunction(uflSpace)
un  = ufl.Coefficient(uflSpace)
vn  = ufl.Coefficient(uflSpace)

# right hand sie (time derivative part + explicit forcing in v)
a_ex = ufl.inner(un, phi) * ufl.dx
# left hand side (heat equation in first variable + backward Euler in time)
a_im = (dt * spiral_D * ufl.inner(ufl.grad(u), ufl.grad(phi)) +
        ufl.inner(u,phi)) * ufl.dx

ustar = (vn[0]+spiral_b)/spiral_a
a_ex += ufl.conditional(un[0]<ustar, dt/spiral_eps* u[0]*(1-un[0])*(un[0]-ustar),
                                     dt/spiral_eps*un[0]*(1-u[0]) *(un[0]-ustar) ) * phi[0] * ufl.dx

equation = a_im == a_ex


# In[ ]:


rhs_gf = create.function("ufl", grid, "rhs", order=2,
                         ufl=ufl.as_vector( [vn[0] + dt*spiral_h(un[0], vn[0]) ]),
                         coefficients={un: uh_n, vn: vh_n} )


# The model is now completely implemented and can be created, together with the corresponding scheme:

# In[ ]:


model = create.model("integrands", grid,  equation,  coefficients={un: uh_n, vn: vh_n} )


# In[ ]:


solverParameters = {
        "fem.solver.newton.tolerance": 1e-3,
        "fem.solver.newton.linabstol": 1e-5,
        "fem.solver.newton.linreduction": 1e-5,
        "fem.solver.newton.verbose": 0,
        "fem.solver.newton.linear.verbose": 0
    }
scheme = create.scheme("galerkin", model, spc, solver="cg",parameters=solverParameters)


# To show the solution we make use of the _animate_ module of _matplotlib_:

# In[ ]:


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
    global t, stepsize, nextstep, iterations, dt
    # print('frame',count,t)
    while t < nextstep:
        uh_n.assign(uh)
        vh_n.assign(vh)
        _,info = scheme.solve(target=uh)
        vh.interpolate( rhs_gf )
        # print("Computing solution a t = " + str(t + dt), "iterations: " + info["linear_iterations"] )
        iterations += int( info["linear_iterations"] )
        t     += dt
    data = uh.pointData(1)
    plt.tricontourf(triangulation, data[:,0], cmap=plt.cm.rainbow, levels=levels)
    grid.writeVTK("spiral", pointdata=[uh], number=count)
    nextstep += stepsize

animation.FuncAnimation(fig, animate, frames=25, interval=100, blit=False)
