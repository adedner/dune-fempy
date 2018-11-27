import numpy
from dune.fem.plotting import plotPointData as plot
import dune.create as create
from dune.grid import structuredGrid
from dune.istl import blockVector

g = structuredGrid([0,0],[1,1],[2,3])

s = create.space("lagrange",g,dimrange=2,storage="istl")
f1 = s.interpolate([2,1], name="tmp")
dofs = blockVector(int(s.size/s.localBlockSize), s.localBlockSize)
# f2 = s.function("tmp", expr=[2,1], dofVector=dofs)
f2 = s.function("tmp", dofVector=dofs)
f2.interpolate([2,1])
assert all([(d1-d2).two_norm==0 for d1,d2 in zip(dofs,f1.as_istl)])
assert all([(d1-d2).two_norm==0 for d1,d2 in zip(dofs,f2.as_istl)])

s = create.space("lagrange",g,dimrange=2,storage="fem")
f1 = s.interpolate([2,1], name="tmp")
dofs = numpy.ndarray(s.size)
f2 = s.function("tmp", [2,1], dofs)
assert not (dofs-f1.as_numpy).any()
assert not (dofs-f2.as_numpy).any()

f1.interpolate([3,2])
dofs[:] = f1.as_numpy
assert not (dofs-f1.as_numpy).any()
assert not (dofs-f2.as_numpy).any()

try:
    import petsc4py
    from petsc4py import PETSc
    s = create.space("lagrange",g,dimrange=2,storage="petsc")
    f1 = s.interpolate([2,1], name="tmp")
    dofs = PETSc.Vec(s.size)
    f2 = s.function("tmp", [2,1], dofs)
    assert all([(d1-d2).two_norm==0 for d1,d2 in zip(dofs,f1.as_petsc)])
    assert all([(d1-d2).two_norm==0 for d1,d2 in zip(dofs,f2.as_petsc)])
except:
    pass
