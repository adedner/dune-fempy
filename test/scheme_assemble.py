import time,sys
import numpy as np
import scipy.sparse
import dune.grid
import dune.fem
import dune.create as create

from dune.ufl import Space
from ufl import TestFunction, TrialFunction, SpatialCoordinate, ds, dx, inner, grad

test_fem   = True
test_eigen = False # not working at the moment
test_istl  = True
test_petsc = True
testLoop = 10
# testLoop = 1000

# grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [89, 89]), dimgrid=2)
grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [9, 9]), dimgrid=2)

d = 0.001
p = 1.7

uflSpace = Space((2,2),1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

rhs = (x[0] + x[1]) * v[0]
a = (pow(d + inner(grad(u), grad(u)), (p-2)/2)*inner(grad(u), grad(v)) + grad(u[0])[0]*v[0]) * dx + 10*inner(u, v) * ds
b = rhs * dx + 10*rhs * ds
model = create.model("integrands", grid, a==b)

def test(space):
    if test_fem:
        femSpace  = create.space(space, grid, dimrange=1, order=1, storage='fem')
        femScheme = create.scheme("galerkin", model, femSpace)
        fem_h     = create.function("discrete", femSpace, name="fem")
        fem_dofs  = fem_h.as_numpy
        # femScheme.solve(target = fem_h)
        start= time.clock()
        for i in range(testLoop):
            fem_mat = femScheme.assemble(fem_h).as_numpy
        end= time.clock()
        # print( "fem:", (end-start)/testLoop, flush=True )
        sys.stdout.flush()
        fem_coo   = fem_mat.tocoo()
        # for i,j,v in zip(fem_coo.row,fem_coo.col,fem_coo.data):
        #     print(i,j,v)
        # print("****************************",flush=True)

    if test_eigen:
        eigenSpace  = create.space(space, grid, dimrange=1, order=1, storage='eigen')
        eigenScheme = create.scheme("galerkin", model, eigenSpace)
        eigen_h     = create.function("discrete", eigenSpace, name="eigen")
        eigen_dofs  = eigen_h.as_numpy
        # eigenScheme.solve(target = eigen_h)
        start= time.clock()
        for i in range(testLoop):
           eigen_mat = eigenScheme.assemble(eigen_h).as_numpy
        end= time.clock()
        # print( "eigen:", (end-start)/testLoop, flush=True )
        sys.stdout.flush()
        eigen_coo = eigen_mat.tocoo()
        # for i,j,v in zip(eigen_coo.row,eigen_coo.col,eigen_coo.data):
        #     print(i,j,v)
        # print("****************************",flush=True)

    if test_istl:
        istlSpace  = create.space(space, grid, dimrange=1, order=1, storage='istl')
        istlScheme = create.scheme("galerkin", model, istlSpace)
        istl_h     = create.function("discrete", istlSpace, name="istl")
        istl_dofs  = istl_h.as_istl
        # istlScheme.solve(target = istl_h)
        start= time.clock()
        for i in range(testLoop):
           istl_mat = istlScheme.assemble(istl_h).as_istl
        end= time.clock()
        # print( "istl:", (end-start)/testLoop, flush=True )
        sys.stdout.flush()
        # there is no way yet to go from istl to scipy - would be nice to have
        # istl_coo = istl_mat.tocoo()
        # for i,j,v in zip(eigen_coo.row,eigen_coo.col,eigen_coo.data):
        #     print(i,j,v)
        # print("****************************",flush=True)

    if test_petsc:
        import petsc4py
        from petsc4py import PETSc
        petsc4py.init(sys.argv)
        petscSpace  = create.space(space, grid, dimrange=1, order=1, storage='petsc')
        petscScheme = create.scheme("galerkin", model, petscSpace)
        petsc_h     = create.function("discrete", petscSpace, name="petsc")
        petsc_dofs  = petsc_h.as_petsc
        # petscScheme.solve(target = petsc_h)
        petsc_mat = petscScheme.assemble(petsc_h).as_petsc
        rptr, cind, vals = petsc_mat.getValuesCSR()
        petsc_coo = scipy.sparse.csr_matrix((vals,cind,rptr)).tocoo()
        start= time.clock()
        for i in range(testLoop):
            petsc_mat = petscScheme.assemble(petsc_h).as_petsc
        end= time.clock()
        print( "petsc:", (end-start)/testLoop, flush=True )
        sys.stdout.flush()

        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_WORLD)
        ksp.setType("cg")
        # ksp.getPC().setType("icc")
        petsc_h.clear()
        res = petsc_h.copy()
        petscScheme(petsc_h, res)
        petsc_mat = petscScheme.assemble(petsc_h).as_petsc
        ksp.setOperators(petsc_mat,petsc_mat)
        ksp.setFromOptions()
        ksp.solve(res.as_petsc, petsc_h.as_petsc)


        print("****************************",flush=True)
        # print(petsc_mat.size, petsc_mat.getSize(), petsc_mat.getSizes())
        # print(petsc_mat.getType())
        # print(type(petsc_mat))
        # print(petscSpace.size)
        # print(petsc_mat.assembled)
        # rptr, cind, vals = petsc_mat.getValuesCSR()
        # petsc_coo = scipy.sparse.csr_matrix((vals,cind,rptr),shape=(100,100)).tocoo()
        # for i,j,v in zip(petsc_coo.row,petsc_coo.col,petsc_coo.data):
        #     print(i,j,v)

    try:
        # print(eigen_coo.data-fem_coo.data)
        assert (eigen_coo.row == fem_coo.row).all()
        assert (eigen_coo.col == fem_coo.col).all()
        assert np.allclose( eigen_coo.data, fem_coo.data )
    except:
        pass

    try:
        # print(eigen_coo.data-fem_coo.data)
        assert (petsc_coo.row == fem_coo.row).all()
        assert (petsc_coo.col == fem_coo.col).all()
        assert np.allclose( petsc_coo.data, fem_coo.data )
    except:
        pass

test("lagrange")
test("dglagrange")
# test("dglegendre") # only works with cube grids
