from __future__ import print_function, division

from ufl import as_vector, dot, grad, cos, pi, SpatialCoordinate, triangle
from dune.grid import structuredGrid, gridFunction
from dune.fem.space import lagrange, combined
from dune.fem.function import tupleDiscreteFunction

grid  = structuredGrid([0, 0], [1, 1], [16, 16])
spc = lagrange(grid, dimrange=1, order=1)
test = spc.interpolate([0],name="test")
space = combined( spc,spc )
x = SpatialCoordinate(triangle)
exact = as_vector([cos(2.*pi*x[0])*cos(2.*pi*x[1]),dot(x,x)])
solution = space.interpolate(exact, name="solution")
df = tupleDiscreteFunction(spc, spc, components=["p", "s"])
#bug assert df[0] == df.p
#bug assert df[1] == df.s
assert len(df.components[0].dofVector) == len(test.dofVector)
df.interpolate(solution)
df.components[0].interpolate(solution[0])
df.p.interpolate(solution[0])

spc3 = lagrange(grid, dimrange=3, order=1)
df = tupleDiscreteFunction(space,spc3,components=["x","y"])
df.x.interpolate(solution)
df.y.interpolate([0,0,0])
df.components[1].interpolate([1,1,1])
#@gridFunction(grid)
# def f(x):
#     r = x*x
#     return [r,r**2,r**3,r**4]
# df.interpolate(f)
df.interpolate([0,1,2,3,4])
