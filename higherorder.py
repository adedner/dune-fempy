
# coding: utf-8

# In[ ]:


from math import pi, log, sqrt
import matplotlib.pyplot as pyplot
import dune.fem
from dune.fem.plotting import plotPointData as plot
from dune.fem.function import integrate
import dune.create as create
grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [4, 4]), dimgrid=2)

from dune.ufl import Space
from ufl import TestFunction, TrialFunction, SpatialCoordinate, cos, as_vector, dx, grad, inner
uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
f = (8*pi*pi+1)*cos(2*pi*x[0])*cos(2*pi*x[1])
exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )
equation = (inner(grad(u), grad(v)) + inner(u,v)) * dx == f * v[0] * dx

fig = pyplot.figure(figsize=(12,12))
for p in range(1,5):
    space = create.space("lagrange", grid, dimrange=1, order=p)
    scheme = create.scheme("h1", space, equation)
    uh,_ = scheme.solve()
    eh = uh - exact
    error = sqrt( integrate(grid, inner(eh, eh) + inner(grad(eh), grad(eh)) , 5)[0] )
    print("order:",p,error)
    plot(uh, figure=(fig, 240+p), colorbar=False, gridLines="black", level=p)
pyplot.show(),
pyplot.close('all')

