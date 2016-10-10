from __future__ import print_function

from mpi4py import MPI

import math
import ufl

import dune.ufl

from dune.fem.view import geometryGridView

import dune.create as create

order=3

# Basic setup
# -----------
# set up reference domain
grid = create.grid("ALUSimplex", "../data/sphere.dgf", dimgrid=2, dimworld=3)
# grid.hierarchicalGrid.globalRefine(1)

# space for discrete solution on Gamma(t)
spc       = create.space("Lagrange", grid, dimrange=grid.dimWorld, order=order)
# define discrete function for Gamma(t) and setup surface grid
positions = spc.interpolate(\
#               lambda x: x*2.,
              lambda x: x * (1.+0.5*math.sin(2.*math.pi*x[0]*x[1])*math.cos(math.pi*x[2])),
              name="positions")

# space for discrete solution on Gamma(t)
surface   = geometryGridView(positions)
spc       = create.space("Lagrange", surface, dimrange=surface.dimWorld, order=order)
# final time and time step
endTime   = 0.25
dt        = 0.0025
theta     = 0.5

# now set up schemes for left and right hand side
# -----------------------------------------------
# u^{n+1} and forcing
solution  = spc.interpolate(lambda x: x, name="solution")
forcing   = spc.interpolate([0,]*surface.dimWorld, name="solution")

# set up left and right hand side models
# --------------------------------------
uflSpace = dune.ufl.Space((surface.dimGrid, surface.dimWorld), surface.dimWorld)
u = ufl.TrialFunction(uflSpace)
v = ufl.TestFunction(uflSpace)
u_n = dune.ufl.GridCoefficient(solution)

a_im = (dt * theta * ufl.inner(ufl.grad(u), ufl.grad(v)) + ufl.inner(u, v)) * ufl.dx
a_ex = (-dt * (1-theta) * ufl.inner(ufl.grad(u), ufl.grad(v)) + ufl.inner(u, v)) * ufl.dx
lhsModel = create.model("elliptic", surface, a_im == 0)
rhsModel = create.model("elliptic", surface, a_ex == 0)

# left hand side scheme
solver    = create.scheme("h1", solution, lhsModel, "lhs")
rhs       = create.scheme("h1", solution, rhsModel, "rhs")

# time loop
# ---------
count   = 0
t       = 0.
surface.writeVTK("mcf"+str(order)+"-0-", pointdata=[solution], number=count)
surface.writeVTK("mcf"+str(order)+"-", pointdata=[solution], number=count, subsampling=3)

def calcRadius(surface):
  # compute R = int_x |x| / int_x 1
  R   = 0
  vol = 0
  for e in surface.elements():
      rule = grid._module.quadratureRule(e.type, 4)
      for p in rule:
          geo = e.geometry
          R   += geo.position(p.position).two_norm * geo.volume * p.weight
          vol += geo.volume * p.weight
  return R/vol

R0 = calcRadius(surface)

while t < endTime:
    rhs(solution, forcing)
    solver.solve(forcing, solution )
    t     += dt
    count += 1
    print("time: ", t)
    surface.writeVTK("mcf"+str(order)+"-0-", pointdata=[solution], number=count)
    surface.writeVTK("mcf"+str(order)+"-3-", pointdata=[solution], number=count, subsampling=3)
    positions.assign(solution.dofVector())
    if count % 1:
        R      = calcRadius( surface )
        Rexact = math.sqrt(R0*R0-4.*t)
        print("R_h=",R, "Rexact=",Rexact, "difference=",R-Rexact)
