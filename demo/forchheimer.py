# implementation of forchheimer equation from https://arxiv.org/pdf/1508.00294.pdf

import dune.fem
from dune.fem.plotting import plotPointData as plot
import dune.create as create
from math import log
from ufl import TestFunction, TrialFunction, SpatialCoordinate, Coefficient, Constant, triangle
from ufl import as_vector, dx, grad, inner, replace, exp, sqrt
from dune.ufl import Space

dune.fem.parameter.append('../data/parameter')
grid = create.grid('ALUConform', dune.grid.cartesianDomain([0, 0], [1, 1], [4, 4]))
space = create.space('lagrange', grid, dimrange=1, order=2, storage='istl')

uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
rho = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
rho_n = Coefficient(uflSpace)   # previous rho
dt = Constant(triangle)         # time step
t  = Constant(triangle)         # current time

abs_drho = sqrt(inner(grad(rho), grad(rho)))
K = 2/(1 + sqrt(1 + 4*abs_drho))
a = (inner((rho - rho_n)/dt, v) + inner(K*grad(rho), grad(v)))*dx
initial = 1/2*(x[0]*x[0] + x[1]*x[1]) - 1/3*(pow(x[0],3) - pow(x[1],3)) + 1
exact = as_vector( [exp(-2*t)*(initial - 1) + 1] )
b = replace(a, {rho: exact})

initial_gf = create.function('ufl', grid, 'initial', 5, initial)
rho_h = space.interpolate(initial_gf, name='rho_h')
rho_h_n = rho_h.copy()

model = create.model('elliptic', grid, a == b, coefficients={rho_n: rho_h_n})

exact_end = as_vector( [exp(-2)*(initial - 1) + 1] )
exact_gf = create.function('ufl', grid, 'exact', 5, exact_end)
def l2error(en, x):
    val = rho_h.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
    return [ val[0]*val[0] ];
l2error_gf = create.function('local', grid, 'error', 5, l2error)

scheme = create.scheme('h1', space, model)

timeStep = 0.05
model.setConstant(dt, timeStep)
error = 0
for eocLoop in range(5):
    print('# step:', eocLoop, ', size:', grid.size(0))
    time = 0.0
    endTime = 1.0
    while time < (endTime + 1e-6):
        model.setConstant(t, time)
        rho_h_n.assign(rho_h)
        scheme.solve(target=rho_h)
        time += timeStep
    error_old = error
    error = sqrt(l2error_gf.integrate()[0])
    if eocLoop == 0:
        eoc = '-'
    else:
        eoc = log(error/error_old)/log(0.5)
    print('|rho_h - rho| =', error, ', eoc =', eoc)
    grid.writeVTK('forchheimer', pointdata=[rho_h, l2error_gf], number=eocLoop)
    grid.hierarchicalGrid.globalRefine(2)

plot(rho_h)
exit()
