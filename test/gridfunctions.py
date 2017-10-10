import math
import numpy
import ufl

import dune.ufl
import dune.create

from dune.grid import cartesianDomain
grid = dune.create.grid("Yasp", cartesianDomain([0,0],[1,1],[4,4]))

# define global grid function

def expr_global(x):
    return [math.sin(x[0])**2, math.sin(x[0])*math.cos(10*x[1]), math.cos(10*x[1])**2]
gf_global = dune.create.function("global", grid, "expr_global", 3, expr_global)

# define local grid function

def expr_local(e,x):
    y = e.geometry.position(x)
    return [gf_global.localFunction(e).evaluate(x)[1] - y.two_norm,10]
gf_local = dune.create.function("local", grid, "expr_local", 3, expr_local)

# interpolate first component of gf_global into Lagrange space

scalar_space = dune.create.space("Lagrange", grid, dimrange=1, order=2)
df_global0 = scalar_space.interpolate(gf_global[0], name="gf_global")

# create function from UFL expression

uflSpace = dune.ufl.Space(grid, 2, field="double")
x = ufl.SpatialCoordinate(ufl.triangle)
const = ufl.Constant(ufl.triangle)
gf_ufl = dune.create.function("ufl", grid, "ufl", 1, ufl.as_vector([ufl.cos(const*x[1]), ufl.sin(x[0])]))

# interpolate ufl function

gf_ufl.setConstant(const, [10])
df_ufl1 = scalar_space.interpolate(gf_ufl[1], name="gf_ufl")

# create function from UFL expression containing coefficient

scalarUFLSpace = dune.ufl.Space(grid,1)
coeffFunc = dune.create.function("global", grid, "global_velocity", 0, lambda x: [2])
gf_ufl = dune.create.function("ufl", grid, "ufl", 4, ufl.as_vector([ufl.sin(x[0])**2 * gf_global[0]]))

# interpolate ufl function into numpy array

scalar_space = dune.create.space("Lagrange", grid, dimrange=1, order=1, storage="eigen")
df_ufl = scalar_space.interpolate(gf_ufl)
#dofs = df_ufl.as_numpy
