# implementation of forchheimer equation from https://arxiv.org/pdf/1508.00294.pdf

import dune.fem
from dune.fem.plotting import plotPointData as plot
import dune.create as create
from math import log
from ufl import TestFunction, TrialFunction, SpatialCoordinate, Coefficient, Constant, triangle
from ufl import as_vector, dx, grad, inner, replace, exp, sqrt
from dune.ufl import Space

dune.fem.parameter.append("../data/parameter")
grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [4, 4]), dimgrid=2)
space = create.space("lagrange", grid, dimrange=1, order=2, storage='istl')

uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
un = Coefficient(uflSpace)   # coefficient used in ufl for old uh
dt = Constant(triangle)      # set to time step size later on
t  = Constant(triangle)      # current time

abs_dphi = sqrt(inner(grad(u), grad(u)))
K = 2/(1 + sqrt(1 + 4*abs_dphi))
a = (inner((u - un)/dt, v) + inner(K*grad(u), grad(v)))*dx
initial = 0.5*(x[0]*x[0] + x[1]*x[1]) - 1/3*(pow(x[0], 3) - pow(x[1], 3)) + 1
exact = as_vector( [exp(-2*t)*(initial - 1) + 1] )
exact_end = as_vector( [exp(-2)*(initial - 1) + 1] )
b = replace(a, {u: exact})

initial_gf = create.function("ufl", grid, "initial", 5, initial)
exact_gf = create.function("ufl", grid, "exact", 5, exact_end)
uh = space.interpolate(initial_gf, name="uh")
uh_n = uh.copy()

def l2error(en, x):
    val = uh.localFunction(en).evaluate(x) - exact_gf.localFunction(en).evaluate(x)
    return [ val[0]*val[0] ];
l2error_gf = create.function("local", grid, "error", 5, l2error)

model = create.model("elliptic", grid, a == b, coefficients={un: uh_n})

scheme = create.scheme("h1", space, model)

timeStep = 0.05
model.setConstant(dt, timeStep)
time = 0.0
endTime = 1.0

error = 10
for eocLoop in range(3):
    print('# step:', eocLoop, ", size:", grid.size(0))
    time = 0.0
    while time < (endTime + 1e-6):
        print(time, grid.size(0), end="\r")
        model.setConstant(t, time)
        uh_n.assign(uh)
        scheme.solve(target=uh)
        time += timeStep
    error_old = error
    error = sqrt(l2error_gf.integrate()[0])
    if eocLoop == 0:
        eoc = 'n/a'
    else:
        eoc = log(error/error_old)/log(0.5)
    plot(uh)
    print("|uh - u| =", error, ", eoc =", eoc)
    #grid.writeVTK('forchheimer', pointdata=[uh, error], number=eocLoop)
    grid.hierarchicalGrid.globalRefine(1)

exit()
