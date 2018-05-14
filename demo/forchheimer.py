# implementation of forchheimer equation from https://arxiv.org/pdf/1508.00294.pdf

import dune.fem
from dune.fem.plotting import plotPointData as plot
import dune.create as create
from math import log
from ufl import TestFunction, TrialFunction, SpatialCoordinate, Coefficient, Constant, triangle
from ufl import as_vector, dx, grad, inner, replace, exp, sqrt, dot
from dune.ufl import Space

grid = create.grid('ALUConform', dune.grid.cartesianDomain([0, 0], [1, 1], [4, 4]))
space = create.space('lagrange', grid, dimrange=1, order=2, storage='istl')

ufl_space = Space((grid.dimGrid, grid.dimWorld), 1)
rho = TrialFunction(ufl_space)
v = TestFunction(ufl_space)
x = SpatialCoordinate(ufl_space.cell())
dt = Constant(triangle)         # time step
t  = Constant(triangle)         # current time

initial = 1/2*(x[0]**2 + x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1
rho_h = space.interpolate(initial, name='rho_h')
rho_h_n = rho_h.copy()

exact_end = as_vector( [exp(-2)*(initial - 1) + 1] )
l2error_gf = create.function('ufl', grid, 'error', 5, dot(rho_h-exact_end, rho_h-exact_end))
h1error_gf = create.function('ufl', grid, 'error', 5, \
                            inner(grad(rho_h-exact_end), grad(rho_h-exact_end)))

abs_drho = sqrt(inner(grad(rho), grad(rho)))
K = 2/(1 + sqrt(1 + 4*abs_drho))
a = (inner((rho - rho_h_n)/dt, v) + inner(K*grad(rho), grad(v)))*dx
exact = as_vector( [exp(-2*t)*(initial - 1) + 1] )
b = replace(a, {rho: exact})

model = create.model('elliptic', grid, a == b)
scheme = create.scheme('h1', space, model)

timeStep = 0.05
model.setConstant(dt, timeStep)
l2error = 0
h1error = 0
for eocLoop in range(5):
    print('# step:', eocLoop, ', size:', grid.size(0))
    time = 0.0
    endTime = 1.0
    rho_h.interpolate(initial)
    while time < (endTime + 1e-6):
        model.setConstant(t, time)
        rho_h_n.assign(rho_h)
        scheme.solve(target=rho_h)
        time += timeStep
    l2error_old = l2error
    h1error_old = h1error
    l2error = sqrt(l2error_gf.integrate()[0])
    h1error = sqrt(h1error_gf.integrate()[0])
    if eocLoop == 0:
        l2eoc = '-'
        h1eoc = '-'
    else:
        l2eoc = log(l2error/l2error_old)/log(0.5)
        h1eoc = log(h1error/h1error_old)/log(0.5)
    print('|rho_h - rho| =', l2error, ', eoc =', l2eoc, \
         ', |grad(rho_h - rho)| =', h1error, ', eoc =', h1eoc)
    plot(rho_h)
    grid.writeVTK('forchheimer', pointdata=[rho_h, l2error_gf, h1error_gf], number=eocLoop)
    grid.hierarchicalGrid.globalRefine(2)

exit()
