# coding: utf-8

# # Solving an Elliptic PDE
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

# In[ ]:

import dune.fem
dune.fem.parameter.append("parameter")


# First, we create our computational grid. Our domain will be the unit square divided into 8x8 quadrilaterals. To actually create the grid, we choose an implementation of the DUNE grid interface: a 2-dimensional ALUGrid with simplices and conforming bisection refinement.

# In[ ]:

import dune.create as create
grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8]), dimgrid=2)


# We set up the base variables u, v and x in UFL.

# In[ ]:

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
# We take $f = 9\pi^2\cos(2\pi x_0)\cos(2\pi x_1)$.
#
# Note that then the exact solution is then
# $u = \cos(2\pi x_0)\cos(2\pi x_1)$.

# In[ ]:

from math import pi
from ufl import cos, as_vector, dx, grad, inner
f = 9*pi*pi*cos(2*pi*x[0])*cos(2*pi*x[1])
exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )
equation = (inner(grad(u), grad(v)) + inner(u,v)) * dx == f * v[0] * dx
equation


# We create the space and the model.

# In[ ]:

spc = create.space("Lagrange", grid, dimrange=1, order=1)
model = create.model("elliptic", grid, equation)


# We create the scheme and set parameters for the solver.

# In[ ]:

scheme = create.scheme("h1", spc, model,       parameters=       {"fem.solver.newton.linabstol": 1e-10,
        "fem.solver.newton.linreduction": 1e-10,
        "fem.solver.newton.verbose": 0,
        "fem.solver.newton.linear.verbose": 0})


# We create a grid function for our exact solution.

# In[ ]:

exact_gf = create.function("ufl", grid, "exact", 5, exact)


# We set up a function for plotting the data using matplotlib.

# In[ ]:

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


# Now we solve the system. We assign the solution to `uh`, and define a function to calculate the $L^2$ error, i.e. $|u_h - u|_{L^2}$. We output the data to a vtk file with name `laplace`, and plot it using `plot`. Finally we refine the grid twice and repeat the process.

# In[ ]:

from math import sqrt
for i in range(2):
    print("solve on level", i, "number of dofs=", grid.size(2))
    uh,_ = scheme.solve()
    def l2error(en,x):
        val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
        return [ val[0]*val[0] ];
    l2error_gf = create.function("local", grid, "error", 5, l2error)
    error = sqrt(l2error_gf.integrate()[0])

    print("size:", grid.size(0), "L2-error:", error)
    grid.writeVTK("laplace", pointdata=[uh, l2error_gf])

    plot(grid, uh)

    if i < 1:
        grid.hierarchicalGrid.globalRefine(2)


# Congratulations! You have successfully solved and visualized your first PDE using dune-fempy.
