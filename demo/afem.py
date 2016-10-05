from __future__ import print_function
import math
from mpi4py import MPI

from ufl import *
from dune.ufl import Space as UFLSpace

import dune.create as create
import dune.grid as grid
import dune.fem as fem
from dune.fem.view import adaptiveLeafGridView


# set the angle for the corner (0<angle<=360)
cornerAngle = 360.
order = 2

# exact solution for this angle
def exact(x):
    r2 = x.two_norm2
    phi = math.atan2(x[1], x[0])
    if x[1] < 0:
        phi += 2*math.pi
    return [(r2**(90./cornerAngle)) * sin(180./cornerAngle*phi)]
def exactJac(x):
    r2 = x.two_norm2
    phi = math.atan2(x[1], x[0])
    if x[1] < 0:
        phi += 2*math.pi
    r2dx=2.*x[0]
    r2dy=2.*x[1]
    if r2 == 0:
        phidx = 0
        phidy = 0
        r2pow = 0
    else:
        phidx=-x[1]/r2
        phidy=x[0]/r2
        r2pow = r2**(90./cornerAngle-1)
    dx = r2pow * ( 90./cornerAngle*r2dx * sin(180./cornerAngle * phi)
                  + r2 * 180./cornerAngle*cos( 180./cornerAngle * phi) * phidx )
    dy = r2pow * ( 90./cornerAngle*r2dy * sin(180./cornerAngle * phi)
                  + r2 * 180./cornerAngle*cos( 180./cornerAngle * phi) * phidy )
    return [dx,dy]

# define the grid for this domain (vertices are the origin and 4 equally spaces on the
# unit sphere starting with (1,0) and ending at # (cos(cornerAngle),sin(cornerAngle))
vertices = [(0,0)]
dgf = "VERTEX\n 0 0 \n"
for i in range(0,5):
    dgf += str(math.cos(cornerAngle/4*math.pi/180*i)) + " " + str(math.sin(cornerAngle/4*math.pi/180*i)) + "\n"
dgf += "#\n"
dgf += """
SIMPLEX
  0 1 2
  0 2 3
  0 3 4
  0 4 5
#
BOUNDARYDOMAIN
  default 1
#
PROJECTION
  function p(x) = x / |x|
  segment  1 2  p
  segment  2 3  p
  segment  3 4  p
  segment  4 5  p
#
"""
grid = create.grid("ALUConform", grid.string2dgf(dgf), dimgrid=2)
grid = adaptiveLeafGridView(grid)
grid.hierarchicalGrid.globalRefine(2)
exact_gf = create.function("global", grid, "exact", order+1, exact)

# use a piecewise quadratic Lagrange space
spc  = create.space( "Lagrange", grid, dimrange=1, order=order)

# the model is -laplace u = 0 with Dirichlet boundary conditions
uflSpace = UFLSpace(2, 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
bnd_u = Coefficient(uflSpace)
a = inner(grad(u), grad(v)) * dx

model = create.model("elliptic", grid, a == 0, dirichlet={1:[bnd_u]},
        tempVars=False, coefficients={bnd_u: exact_gf})

# set up the scheme
laplace = create.scheme("h1", spc, model, "afem")
uh = spc.interpolate(lambda x: [0], name="solution")
laplace.solve(target=uh)

# function used for computing approximation error
def h1error(en,x):
    y = en.geometry.position(x)
    val = uh.localFunction(en).evaluate(x) - exact(y)
    jac = uh.localFunction(en).jacobian(x)[0] - exactJac(y)
    return [ sqrt( val[0]*val[0] + jac*jac) ];
h1error_gf = create.function( "local", grid, "error", order+1, h1error )

# adaptive loop (mark, estimate, solve)
count = 0
tol = 0.05 # use 0 for global refinement
while count < 20:
    error = math.sqrt(h1error_gf.integrate()[0])
    [estimate, marked] = laplace.mark(uh, tol)
    grid.writeVTK("afem", pointdata=[uh],
            celldata=[create.function("levels",grid)], number=count )
    print(count, ": size=",grid.size(0), "estimate=",estimate,"error=",error)
    if marked == False or estimate < tol:
        break
    if tol == 0.:
        grid.hierarchicalGrid.globalRefine(2)
        uh.interpolate([0])  # initial guess needed
    else:
        fem.adapt(grid.hierarchicalGrid, [uh])
        fem.loadBalance(grid.hierarchicalGrid, [uh])
    laplace.solve( target=uh )
    count += 1
