import time
import numpy
import math
import dune.plotting
dune.plotting.block = False

from dune.grid import structuredGrid
grid = structuredGrid([0, 0], [1, 1], [4, 4])

import dune.create as create
space = create.space('lagrange', grid, dimrange=1, order=2)

from ufl import SpatialCoordinate
x = SpatialCoordinate(space)

initial = 1/2*(x[0]**2 + x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1
u_h = space.interpolate(initial, name='u_h')
u_h_n = u_h.copy(name="previous")

from ufl import TestFunction, TrialFunction
from dune.ufl import NamedConstant
u = TrialFunction(space)
v = TestFunction(space)
dt = NamedConstant(space, "dt")    # time step
t  = NamedConstant(space, "t")     # current time

from ufl import dx, grad, div, inner, sqrt
abs_du = lambda u: sqrt(inner(grad(u), grad(u)))
K = lambda u: 2/(1 + sqrt(1 + 4*abs_du(u)))
a = ( inner((u - u_h_n)/dt, v) \
    + 0.5*inner(K(u)*grad(u), grad(v)) \
    + 0.5*inner(K(u_h_n)*grad(u_h_n), grad(v)) ) * dx

from ufl import as_vector, exp
exact = lambda t: as_vector([exp(-2*t)*(initial - 1) + 1])

from ufl import dot, FacetNormal, ds
n = FacetNormal(space)
b = inner(-2*exp(-2*t)*(initial - 1) \
    - div(K(exact(t))*grad(exact(t)[0])), v[0]) * dx \
    + K(exact(t))*dot(grad(exact(t)[0]), n) * v[0] * ds

scheme = create.scheme("galerkin", a == b, solver='cg')

scheme.model.dt = 0.001

def evolve(scheme, u_h, u_h_n):
    time = 0
    endTime = 1.0
    while time < (endTime - 1e-6):
        scheme.model.t = time + 0.5*scheme.model.dt
        u_h_n.assign(u_h)
        scheme.solve(target=u_h)
        time += scheme.model.dt

exact_end = exact(1)
l2error_fn = inner(u_h - exact_end, u_h - exact_end)
h1error_fn = inner(grad(u_h - exact_end), grad(u_h - exact_end))

from math import log
from dune.fem.function import integrate
from dune.fem.plotting import plotPointData as plot
l2error = 0
h1error = 0
for eocLoop in range(5):
    print('# step:', eocLoop, ', size:', grid.size(0), \
           ', dt:', scheme.model.dt)
    start_time = time.time()
    u_h.interpolate(initial)
    evolve(scheme, u_h, u_h_n)
    print('solve time:', time.time() - start_time)
    l2error_old = l2error
    h1error_old = h1error
    l2error = sqrt( integrate(grid, l2error_fn, 5)[0] )
    h1error = sqrt( integrate(grid, h1error_fn, 5)[0] )
    if eocLoop == 0:
        l2eoc = '-'
        h1eoc = '-'
    else:
        l2eoc = log(l2error/l2error_old)/log(0.5)
        h1eoc = log(h1error/h1error_old)/log(0.5)
    print('|u_h - u| =', l2error, ', eoc =', l2eoc,
          ', |nabla(u_h - u)| =', h1error, ', eoc =', h1eoc)
    #plot(u_h)
    grid.writeVTK('forchheimer', pointdata={'u': u_h, 'l2error':
                  l2error_fn, 'h1error': h1error_fn}, number=eocLoop)
    grid.hierarchicalGrid.globalRefine(1)
    scheme.model.dt /= 2
