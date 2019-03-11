# <markdowncell>
# ## Mean Curvature Flow
#
# In this example we compute the mean curvature flow of a surface:
# \begin{align}
#    \frac{\partial}{\partial_t} x &= H(x)  && \text{for } x\in\Gamma(t)
# \end{align}
# Assume we can define a reference surface $\Gamma_0$ such that
# we can write the evolving surface $\Gamma(t)$ in the form
# \begin{gather}
#   \Gamma(t) = X(t,\Gamma_0)
# \end{gather}
# It is now possible to show that the vector valued function $X=X(t,x)$ with $x\in\Gamma_0$ satisfies:
# \begin{gather}
#   \frac{\partial}{\partial_t}X = - H(X)\nu(X)
# \end{gather}
# where $H$ is the mean curvature of $\Gamma_t$ and $\nu$ is its outward pointing normal.
#
# We will solve this using a finite element approach based on the following time discrete approximation:
# \begin{gather}
#   \int_{\Gamma^n} \big( U^{n+1} - {\rm id}\big) \cdot \varphi +
#     \tau \int_{\Gamma^n} \big(
#     \theta\nabla_{\Gamma^n} U^{n+1} + (1-\theta) I \big)
#     \colon\nabla_{\Gamma^n}\varphi
#   =0~.
# \end{gather}
# Here $U^n$ parametrizes $\Gamma(t^{n+1})$ over
# $\Gamma^n:=\Gamma(t^{n})$,
# $I$ is the identity matrix, $\tau$ is the time step and
# $\theta\in[0,1]$ is a discretization parameter.
# <img src="mcf.gif" style="height:228px;">

# <codecell>
from __future__ import print_function
try:
    get_ipython().magic('matplotlib inline # can also use notebook or nbagg')
except:
    pass

import math, time
import pickle

from ufl import *
import dune.ufl
from dune.generator import algorithm
import dune.geometry as geometry
import dune.fem as fem
from dune.fem.plotting import plotPointData as plot
import matplotlib.pyplot as pyplot
from IPython import display

# polynomial order of surface approximation
order = 2

# initial radius
R0 = 2.

# end time
endTime = 0.1


# <markdowncell>
# Main function for calculating the mean curvature flow of a given surface.
# If first argument is `True` the radius of the computed surface is
# computed using an algorithm implemented in C++ otherwise the computation
# is done in Python.
#
# Timings for a number of different grid refinements is dumped to disk


# <codecell>
from dune.fem.view import geometryGridView as geoGridView
from dune.fem.space import lagrange as solutionSpace
from dune.fem.scheme import galerkin as solutionScheme
def calculate(use_cpp, gridView):
    # space on Gamma_0 to describe position of Gamma(t)
    space = solutionSpace(gridView, dimRange=gridView.dimWorld, order=order)
    u = TrialFunction(space)
    v = TestFunction(space)
    x = SpatialCoordinate(space.cell())
    positions = space.interpolate(x, name="position")

    # space for discrete solution on Gamma(t)
    surface = geoGridView(positions)
    space = solutionSpace(surface, dimRange=surface.dimWorld, order=order)
    solution  = space.interpolate(x, name="solution")

    # set up model using theta scheme
    theta = 0.5   # Crank-Nicholson

    I = Identity(3)
    dt = dune.ufl.Constant(0,"dt")

    a = (inner(u - x, v) + dt * inner(theta*grad(u)
        + (1 - theta)*I, grad(v))) * dx

    scheme = solutionScheme(a == 0, space, solver="cg")

    if use_cpp:
        radius = algorithm.load('calcRadius', 'radius.hh', surface)
        file_path = 'cpp_time.p'
    else:
        # compute an averaged radius of the surface
        def radius(surface):
            # compute R = int_x |x| / int_x 1
            R   = 0
            vol = 0
            for e in surface.elements:
                rule = geometry.quadratureRule(e.type, 4)
                for p in rule:
                    geo = e.geometry
                    weight = geo.volume * p.weight
                    R   += geo.toGlobal(p.position).two_norm * weight
                    vol += weight
            return R/vol
        file_path = 'python_time.p'

    scheme.model.dt = 0.02

    import numpy as np
    pyplot.figure()
    pyplot.gca().set_xlim([0, endTime])
    pyplot.gca().set_ylabel("error")
    pyplot.gca().set_xlabel("time")

    numberOfLoops = 3
    times = np.zeros(numberOfLoops)
    errors = np.zeros(numberOfLoops)
    totalIterations = np.zeros(numberOfLoops, np.dtype(np.uint32))
    gridSizes = np.zeros(numberOfLoops, np.dtype(np.uint32))
    for i in range(numberOfLoops):
        positions.interpolate(lambda x: x * (R0/x.two_norm))
        solution.interpolate(lambda x: x)
        t = 0.
        R = radius( surface )
        Rexact = math.sqrt(R0*R0 - 4.*t)
        x = np.array([t])
        y = np.array([R - Rexact])
        iterations = 0
        start = time.time()
        while t < endTime:
            info = scheme.solve(target=solution)
            # move the surface
            positions.dofVector.assign(solution.dofVector)
            # store some information about the solution process
            iterations += int( info["linear_iterations"] )
            t          += scheme.model.dt
            R           = radius( surface )
            Rexact      = math.sqrt(R0*R0-4.*t)
        print("time used:", time.time() - start)
        times[i] = time.time() - start
        errors[i] = abs(R-Rexact)
        totalIterations[i] = iterations
        gridSizes[i] = gridView.size(2)
        if i < numberOfLoops - 1:
            gridView.hierarchicalGrid.globalRefine(1)
            scheme.model.dt /= 2
    eocs = np.log(errors[0:][:numberOfLoops-1] / errors[1:]) / math.log(math.sqrt(2))
    try:
        import pandas as pd
        keys = {'size': gridSizes, 'error': errors, "eoc": np.insert(eocs, 0, None), 'iterations': totalIterations}
        table = pd.DataFrame(keys, index=range(numberOfLoops),columns=['size', 'error', 'eoc', 'iterations'])
        print(table)
    except ImportError:
        print("pandas could not be used to show table with results")
        pass
    pickle.dump([gridSizes, times], open(file_path,'wb'))


# <markdowncell>
# Compute the mean curvature flow evolution of a spherical surface. Compare
# computational time of a pure Python implementation and using a C++
# algorithm to compute the radius of the surface for verifying the
# algorithm.

# <codecell>
# set up reference domain Gamma_0
from dune.alugrid import aluConformGrid as leafGridView
gridView = leafGridView("sphere.dgf", dimgrid=2, dimworld=3)
calculate(True, gridView)
gridView = leafGridView("sphere.dgf", dimgrid=2, dimworld=3)
calculate(False, gridView)
