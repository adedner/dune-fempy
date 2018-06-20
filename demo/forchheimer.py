# implementation of forchheimer equation from https://arxiv.org/pdf/1508.00294.pdf

import dune.fem
import dune.create as create
from dune.fem.function import integrate
from dune.fem.plotting import plotPointData as plot
from dune.grid import structuredGrid
from math import log
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle
from ufl import as_vector, dx, grad, inner, replace, exp, sqrt, dot
from dune.ufl import NamedConstant, Space

grid = structuredGrid([0, 0], [1, 1], [4, 4])
space = create.space('lagrange', grid, dimrange=1, order=2, storage='istl')

u = TrialFunction(space)
v = TestFunction(space)
x = SpatialCoordinate(space.cell())
dt = NamedConstant(triangle, "dt")    # time step
t  = NamedConstant(triangle, "t")     # current time

initial = 1/2*(x[0]**2 + x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1
u_h = space.interpolate(initial, name='u_h')
u_h_n = u_h.copy()

exact_end = as_vector( [exp(-2)*(initial - 1) + 1] )
l2error_fn = dot(u_h - exact_end, u_h - exact_end)
h1error_fn = inner(grad(u_h - exact_end), grad(u_h - exact_end))

abs_du = sqrt(inner(grad(u), grad(u)))
K = 2/(1 + sqrt(1 + 4*abs_du))
a = (inner((u - u_h_n)/dt, v) + inner(K*grad(u), grad(v)))*dx
exact = as_vector( [exp(-2*t)*(initial - 1) + 1] )
b = replace(a, {u: exact})

scheme = create.scheme("h1", space, a == b, solver='cg')

def evolve(scheme, u_h, u_h_n):
    time = 0
    endTime = 1.0
    while time < (endTime + 1e-6):
        scheme.model.t = time
        u_h_n.assign(u_h)
        scheme.solve(target=u_h)
        time += scheme.model.dt

scheme.model.dt = 0.05
l2error = 0
h1error = 0
for eocLoop in range(3):
    print('# step:', eocLoop, ', size:', grid.size(0))
    u_h.interpolate(initial)
    evolve(scheme, u_h, u_h_n)
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
    print('|u_h - u| =', l2error, ', eoc =', l2eoc, \
         ', |grad(u_h - u)| =', h1error, ', eoc =', h1eoc)
    plot(u_h)
    grid.writeVTK('forchheimer', pointdata={'u': u_h, 'l2error': l2error_fn,
                  'h1error': h1error_fn}, number=eocLoop)
    grid.hierarchicalGrid.globalRefine(1)

exit()
