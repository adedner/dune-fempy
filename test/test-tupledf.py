from __future__ import print_function, division

from ufl import as_vector, dot, grad, cos, pi, SpatialCoordinate, triangle
from dune.grid import structuredGrid, gridFunction
from dune.fem.space import lagrange, combined
from dune.fem.function import tupleDiscreteFunction

grid  = structuredGrid([0, 0], [1, 1], [16, 16])
spc = lagrange(grid, dimrange=1, order=1)
space = combined( spc,spc )
x = SpatialCoordinate(triangle)
exact = as_vector([cos(2.*pi*x[0])*cos(2.*pi*x[1]),dot(x,x)])
solution = space.interpolate(exact, name="solution")

test = space.interpolate(grad(solution[0]), name="tmp")

df = tupleDiscreteFunction(spc,spc,components=["p","s"])
# assert df[0] == df.p
# assert df[1] == df.s
# bug assert len(df[0].dofVector) == len(solution.dofVector)
df.interpolate(solution)
# df.p.interpolate(solution[0])
df.component(0).interpolate(solution[0])
# bug df[0].interpolate(solution[0])
# bug df[1].interpolate(solution[1])

spc2 = lagrange(grid, dimrange=2, order=1)
df = tupleDiscreteFunction(space,spc2,components=["x","y"])
# assert df[0] == df.x
# assert df[1] == df.y
# assert len(df.x.p.dofVector) == len(solution.dofVector)
# assert len(df.x.p) == len(solution)
# assert len(df.x == df.y)
# df.x[0].interpolate(solution[0])
# df.x[1].interpolate(solution[1])
# df[0][0].interpolate(solution[0])
# df[0][1].interpolate(solution[1])
# df[1].interpolate(solution)
# @gridFunction(grid)
# def f(x):
#     r = x*x
#     return [r,r**2,r**3,r**4]
# df.interpolate(f)
df.interpolate([0,1,2,3])
