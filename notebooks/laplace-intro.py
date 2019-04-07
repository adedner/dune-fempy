
# coding: utf-8

# # Solving an Elliptic PDE [(Notebook)][1]
#
# [1]: _downloads/laplace-intro.ipynb
#
# This demo introduces basic usage of dune-fempy, using the Poisson equation as an example. Namely,
#
# \begin{align*}
#   - \Delta u + u &= f && \text{in } \Omega \\
#   \nabla u \cdot \textbf{n} &= 0 && \text{on } \Gamma
# \end{align*}
#
#
# If you have compiled DUNE against MPI, we strongly advise you to first initialize MPI from Python.
# At least OpenMPI is known to fail, if initialized only in the dune-fempy library.

# In[1]:


try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass
import dune.fem
from dune.fem.plotting import plotPointData as plot
dune.fem.parameter.append({"fem.verboserank": "0"})
dune.fem.parameter.append("parameter")


# First, we create our computational grid. Our domain will be the unit square divided into 8x8 quadrilaterals. To actually create the grid, we choose an implementation of the DUNE grid interface: a 2-dimensional ALUGrid with simplices and conforming bisection refinement.

# In[2]:


import dune.create as create
grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [4, 4]), dimgrid=2)


# We set up the base variables u, v and x in UFL.

# In[3]:


from dune.ufl import Space
from ufl import TestFunction, TrialFunction, SpatialCoordinate
uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())


# Next we define the equation for the weak form, given by
#
# \begin{equation}
# \int_{\Omega} uv + \nabla u\cdot\nabla v \ dx =  \int_{\Omega} f v \ dx.
# \end{equation}
# We take $f = 8\pi^2\cos(2\pi x_0)\cos(2\pi x_1)$.
#
# Note that then the exact solution is then
# $u = \cos(2\pi x_0)\cos(2\pi x_1)$.

# In[4]:


from math import pi,log10
from ufl import cos, as_vector, dx, grad, inner
f = (8*pi*pi+1)*cos(2*pi*x[0])*cos(2*pi*x[1])
exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )
equation = (inner(grad(u), grad(v)) + inner(u,v)) * dx == f * v[0] * dx


# We create the space and the model.

# In[5]:


spc = create.space("lagrange", grid, dimRange=1, order=1)
model = create.model("integrands", grid, equation)


# We create the scheme and set parameters for the solver.

# In[6]:


scheme = create.scheme("galerkin", model, spc, solver="cg",
        parameters=
       {"newton.tolerance": 1e-5, "newton.verbose": True,
        "newton.linear.tolerance": 1e-8,
        "newton.linear.preconditioning.method": "ilu",
        "newton.linear.preconditioning.iterations": 1, "newton.linear.preconditioning.relaxation": 1.2,
        "newton.linear.verbose": True})


# We create a grid function for our exact solution.

# In[7]:


exact_gf = create.function("ufl", grid, "exact", 5, exact)


# Now we solve the system. We assign the solution to `uh`, and define a function to calculate the $L^2$ error, i.e. $|u_h - u|_{L^2}$. We output the data to a vtk file with name `laplace`, and plot it using `plot`. Finally we refine the grid twice and repeat the process.

# In[8]:


from math import sqrt
levels=2
for i in range(levels):
    print("solve on level", i, "number of dofs=", spc.size)
    uh = spc.interpolate([0],name="uh")
    scheme.solve(target=uh)
    def l2error(en,x):
        val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
        return [ val[0]*val[0] ]
    l2error_gf = create.function("local", grid, "error", 5, l2error)
    error = sqrt(l2error_gf.integrate()[0])

    print("size:", grid.size(0), "L2-error:", error)
    grid.writeVTK("laplace", pointdata=[uh, l2error_gf])
    plot(uh,gridLines="black")

    if i < levels-1:
        grid.hierarchicalGrid.globalRefine(2)


# Congratulations! You have successfully solved and visualized your first PDE using dune-fempy.
#
# But it is still rather coarse... so either we refine the grid a bit further or we use a second order method:

# In[9]:


spc = create.space("lagrange", grid, dimRange=1, order=2)
# create the scheme but change some of the default parameters..
scheme = create.scheme("galerkin", model, spc, solver="cg",
        parameters=\
       {"newton.tolerance": 1e-5, "newton.verbose": "true",
        "newton.linear.absolutetol": 1e-8, "newton.linear.reductiontol": 1e-8,
        "newton.linear.preconditioning.method": "ilu",
        "newton.linear.preconditioning.iterations": 1, "newton.linear.preconditioning.relaxation": 1.2,
        "newton.linear.verbose": "true"})

print("solve with second order", i, "number of dofs=", spc.size)

uh = spc.interpolate([0],name="uh")
info = scheme.solve(target=uh)
def l2error(en,x):
    val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
    return [ val[0]*val[0] ];
l2error_gf = create.function("local", grid, "error", 5, l2error)

error = sqrt(l2error_gf.integrate()[0])
print("size:", grid.size(0), "L2-error:", error)
plot(uh,level=5)


# Thats looks better already.
#
# Plot the log differecence between the solution and the exact solution - also retrieve the information returned by the solver:

# In[10]:


def error(en,x):
    val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
    return [ log10(abs(val[0])) ];
error_gf = create.function("local", grid, "error", 5, error)

plot(error_gf, level=5,clim=[-5,-2])
print("Solver information: ", info)
