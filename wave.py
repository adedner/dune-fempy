import time
import math
import dune.fem as fem
from dune.grid import reader
from dune.alugrid import aluConformGrid as leafGridView
from dune.fem.view import adaptiveLeafGridView as adaptiveGridView
from dune.fem.space import lagrange as solutionSpace
from ufl import TrialFunction, TestFunction, grad, dot, dx
from dune.ufl import Constant, BoxDirichletBC, DirichletBC
from dune.fem.function import levelFunction, partitionFunction

fem.parameter.append({"fem.verboserank":-1})

order = 1
dimDomain = 2
domain = (reader.gmsh, "wave_tank.msh")
gridView  = adaptiveGridView( leafGridView( domain, dimgrid=dimDomain ) )
gridView.hierarchicalGrid.loadBalance()
V = solutionSpace(gridView, order=order, storage="fem")

p    = V.interpolate(0,name="p")
phi  = V.interpolate(0,name="phi")

u    = TrialFunction(V)
v    = TestFunction(V)
p_in = Constant(0.0)
# the following is equivalent t
# x    = SpatialCoordinate(V)
# bc   = DirichletBC(V, p_in, conditional(x[0]<1e-10,1,0))
bc   = BoxDirichletBC(V, p_in, [None,None],[0,None], eps=1e-10)
T = 2. # 5. # 10.
dt = 0.001
t = 0

from dune.fem.operator import galerkin as operator
op = operator([dot(grad(u),grad(v))*dx,bc])

lumping = operator(u*v*dx)
lapPhi = V.interpolate(1,name="e")
lumped = lapPhi.copy()
lumping(lapPhi,lumped)

pVec = p.as_numpy
phiVec = phi.as_numpy
lapPhiVec = lapPhi.as_numpy
lumpedVec = lumped.as_numpy
lumpedVec[:] = 1/lumpedVec

vtk = gridView.sequencedVTK("wave", pointdata=[p,phi],
       celldata=[levelFunction(gridView), partitionFunction(gridView)])
step = 0
start = time.time()
while t <= T:
    step += 1
    phiVec[:] = phiVec[:]-pVec[:]*dt/2
    op(phi,lapPhi)
    pVec[:] = pVec[:]+dt*lapPhiVec[:]*lumpedVec[:]
    op.setConstraints(p)
    phiVec[:] = phiVec[:]-pVec[:]*dt/2
    if step % 10 == 0:
        vtk()
    t += dt
    op.model.c0 = math.sin(2*math.pi*5*t)
timing = time.time()-start
print("\n runtime:", timing)
