# NOTE: there is some issue with failing convergence when using solver=cg -
# it should work...
import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, dot, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from dune.fem import parameter

from stokesclass import Stokes

parameter.append({"fem.verboserank": "0"})

order = 2
grid = create.grid("ALUCube",constructor=cartesianDomain([0,0],[3,1],[30,10]))
spcU = create.space("lagrange", grid, dimrange=grid.dimension, order=order, storage="istl")
spcP = create.space("lagrange", grid, dimrange=1, order=order-1, storage="istl")

mu = 0.1
Re = 100

cell  = spcU.cell()
x     = SpatialCoordinate(cell)
exact_u     = as_vector( [x[1] * (1.-x[1]), 0] )
exact_p     = as_vector( [ (-2*x[0] + 2)/Re ] )
f           = mu*exact_u if mu>0 else None

velocity = spcU.interpolate([0,]*spcU.dimRange, name="velocity")
pressure = spcP.interpolate([0], name="pressure")

def plot(count=None):
    grid.writeVTK("Stokes",
            pointdata={"pressure":pressure, "exact_p":exact_p},
            pointvector={"velocity":velocity, "exact_velo":exact_u},
            number=count
    )

bcs = [DirichletBC(spcU,exact_u,1)]
stokes = Stokes([spcU,spcP],f,bcs,Re)
stokes.prepare(mu=mu)
stokes.solve(target=[velocity,pressure])

plot()
