# <markdowncell>
# # Mean Curvature Flow (revisited)
#
# We simulate the shrinking of a sphere under mcf using a finite element
# approach based on
# the following time discrete approximation:
# \begin{align}
# \int_{\Gamma^n} \big( U^{n+1} - {\rm id}\big) \cdot \varphi +
#     \tau \int_{\Gamma^n} \big(
#     \theta\nabla_{\Gamma^n} U^{n+1} + (1-\theta) I \big)
#     \colon\nabla_{\Gamma^n}\varphi
#   =0~.
# \end{align}
# Here $U^n$ parametrizes $\Gamma(t^{n+1})$ over
# $\Gamma^n:=\Gamma(t^{n})$,
# $I$ is the identity matrix, $\tau$ is the time step and
# $\theta\in[0,1]$ is a discretization parameter.
#
# In this case we have an exact formula for the evolution of the radius of
# the surface
#
# $$R(t) = \sqrt{R(0)^2 - 4t}.$$
#
# To compare the accuracy of the surface approximation we compute an
# average radius of the discrete surface in each time step $t^n$ using
#
# $$R_h^n = \frac{ \int_{\Gamma^n} |x| }{ |\Gamma^n| }.$$
#
# Computing $R_h^n$ requires a grid traversal and a number of calls to
# interface methods on each element. Doing this on the Python side has a
# potential performance impact which we investigate here by comparing a
# pure python implementation with a hybrid approach where computing $R_h^n$
# is implemented in C++ using the `dune.generator.algorithm` functionality.

# <codecell>
import math, time
import pickle
import numpy as np

from ufl import *
import dune.ufl
from dune.generator import algorithm
import dune.geometry as geometry
import dune.fem as fem
from mcf_cmp_plot import plot

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
# <codecell>
def calcRadius(surface):
    R,vol = 0, 0
    for e in surface.elements:
        rule = geometry.quadratureRule(e.type, 4)
        for p in rule:
            geo = e.geometry
            weight = geo.volume * p.weight
            R   += geo.toGlobal(p.position).two_norm * weight
            vol += weight
    return R/vol

switchCalcRadius = lambda use_cpp,surface: \
             algorithm.load('calcRadius', 'radius.hh', surface) \
             if use_cpp else calcRadius

# <markdowncell>
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

    Rexact = lambda t: math.sqrt(R0*R0 - 4.*t)
    radius = switchCalcRadius(use_cpp,surface)

    scheme.model.dt = 0.02

    numberOfLoops = 3
    times = np.zeros(numberOfLoops)
    errors = np.zeros(numberOfLoops)
    totalIterations = np.zeros(numberOfLoops, np.dtype(np.uint32))
    gridSizes = np.zeros(numberOfLoops, np.dtype(np.uint32))
    for i in range(numberOfLoops):
        positions.interpolate(x * (R0/sqrt(dot(x,x))))
        solution.interpolate(x)
        t = 0.
        error = abs(radius(surface)-Rexact(t))
        iterations = 0
        start = time.time()
        while t < endTime:
            info = scheme.solve(target=solution)
            # move the surface
            positions.dofVector.assign(solution.dofVector)
            # store some information about the solution process
            iterations += int( info["linear_iterations"] )
            t          += scheme.model.dt
            error       = max(error, abs(radius(surface)-Rexact(t)))
        print("time used:", time.time() - start)
        times[i] = time.time() - start
        errors[i] = error
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
    return gridSizes, times


# <markdowncell>
# Compute the mean curvature flow evolution of a spherical surface. Compare
# computational time of a pure Python implementation and using a C++
# algorithm to compute the radius of the surface for verifying the
# algorithm.

# <codecell>
# set up reference domain Gamma_0
results = []
from dune.alugrid import aluConformGrid as leafGridView
gridView = leafGridView("sphere.dgf", dimgrid=2, dimworld=3)
results += [calculate(True, gridView)]

gridView = leafGridView("sphere.dgf", dimgrid=2, dimworld=3)
results += [calculate(False, gridView)]

# <markdowncell>
# Compare the hybrid and pure Python versions

# <codecell>
plot(results[0],results[1])
