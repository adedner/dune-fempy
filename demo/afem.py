from __future__ import print_function
import math
from mpi4py import MPI

from ufl import *
from dune.ufl import Space as UFLSpace
from dune.models.elliptic import compileUFL, importModel
import dune.fem as fem

# set the angle for the corner (0<angle<=360)
cornerAngle = 360

# exact solution for this angle
def exact(x):
    phi = math.atan2(x[1], x[0])
    if x[1] < 0:
        phi += 2*math.pi
    return [(x.two_norm2**(90./cornerAngle)) * sin(180./cornerAngle*phi)]

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
grid = fem.leafGrid(fem.string2dgf(dgf), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
grid.globalRefine(2)
exact_gf = grid.globalGridFunction("exact", exact)

# use a piecewise quadratic Lagrange space
spc  = fem.create.space( "Lagrange", grid, dimrange=1, polorder=2)

# the model is -laplace u = 0 with Dirichlet boundary conditions
uflSpace = UFLSpace(2, 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
bnd_u = Coefficient(uflSpace)
a = inner(grad(u), grad(v)) * dx

model = importModel(grid, a == 0, dirichlet={1:[bnd_u]}, tempVars=False).get()
model.setCoefficient(bnd_u, exact_gf)

# set up the scheme
laplace = fem.create.scheme("FemScheme", spc, model, "afem")
uh = spc.interpolate(lambda x: [0])
laplace.solve(target=uh)

# adaptive loop (mark, estimate, solve)
count = 0
tol = 0.05
while count < 20:
    error = grid.distance(uh,exact_gf)
    [estimate, marked] = laplace.mark(uh, tol)
    grid.writeVTK("afem", pointdata=[uh], celldata=[grid.levelFunction()], number=count )
    print(count, ": size=",grid.size(0), "estimate=",estimate,"error=",error)
    if marked == False or estimate < tol:
        break
    grid.hierarchicalGrid.globalRefine(2)
    uh.interpolate([0])
    # grid.hierarchicalGrid.adapt([uh])
    # grid.hierarchicalGrid.loadBalance([uh])
    laplace.solve( target=uh )
    count += 1
