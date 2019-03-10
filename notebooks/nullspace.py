import math,sys
try:
    import petsc4py, sys
    from petsc4py import PETSc
    petsc4py.init(sys.argv)
except:
    print("Example needs petsc4py")
    sys.exit(0)

import ufl
import dune.create as create
import dune.fem as fem
import dune.geometry
from dune.fem import parameter
from dune.fem.function import integrate

order = 1
surface = create.grid("ALUConform", "sphere.dgf", dimgrid=2, dimworld=3)
spc = create.space("Lagrange", surface, dimRange=1, order=order, storage="petsc")

from ufl import as_vector, TestFunction, TrialFunction, SpatialCoordinate, dx, inner, grad
from dune.ufl import Constant
u = TrialFunction(spc)
v = TestFunction(spc)
x = SpatialCoordinate(spc)
f = as_vector( [(4*x[0]*x[0])-2*((x[2]*x[2])+(x[1]*x[1]))] )
discreteMass = Constant(0,"discreteMass")

equation = inner(grad(u), grad(v)) * dx == (f[0]-discreteMass)*v[0] * dx
scheme = create.scheme("galerkin", equation, spc,
                        parameters={"petsc.preconditioning.method":"hypre"})
exact = as_vector([x[0]*x[0]-1/3])

def correctMass(nullsp, vec):
    Id   = spc.interpolate([1],name="id")
    Area = Id.integrate()[0]
    w    = spc.function("w",dofVector=vec)
    wInt = w.integrate()[0] / Area
    w   -= spc.interpolate([wInt],name="tmp")

def compute():
    Id   = spc.interpolate([1],name="id")
    Area = Id.integrate()[0]
    IntF = integrate(surface,f,order=5)[0]
    scheme.model.discreteMass = IntF/Area

    uh = spc.interpolate([0],name="solution")

    # get matrix (assembled)
    linOp = dune.fem.operator.linear(scheme)
    matrix = linOp.as_petsc

    # Attach nullspace for matrix
    nullsp = PETSc.NullSpace().create(constant=True)
    # nullsp.setFunction(correctMass)     # doesn't work
    matrix.setNullSpace(nullsp)
    matrix.setTransposeNullSpace(nullsp)

    # set up linear solver (using cg and hypre amg preconditioner)
    ksp = PETSc.KSP()
    ksp.create()
    ksp.setType("cg")
    ksp.getPC().setType('hypre') # works well if you included it in your petsc setup
    ksp.setTolerances(rtol=1e-50, atol=1e-15, divtol=None, max_it=10000)
    ksp.setMonitor(lambda ksp,i,r: print(i,r,flush=True))

    # attach matrix to linear solver
    ksp.setOperators(matrix)
    ksp.setFromOptions()

    # compute rhs and project out nullspace
    rhs = uh.copy()
    res = uh.copy()
    rhs.clear()
    scheme(uh, res)
    rhs -= res
    print("Mass of rhs:", rhs.scalarProductDofs(Id))
    # projection of rhs for safty (but hardly changes anything)
    nullsp.remove(rhs.as_petsc)
    print("Mass of rhs after projection:", rhs.scalarProductDofs(Id))

    # solve system
    print("initial residual:", math.sqrt( rhs.as_petsc.dot(rhs.as_petsc) ) )
    ksp.solve(rhs.as_petsc, uh.as_petsc)
    scheme(uh, res)
    print("final residual:", math.sqrt( res.as_petsc.dot(res.as_petsc) ) )

    print("Mass of solution:", uh.integrate()[0],
          " and sum of dofs:", uh.scalarProductDofs(Id))

    # remove mass from final solution
    intUh = uh.integrate()[0]
    uh = spc.interpolate( [uh[0]-intUh/Area], name="solution")
    scheme(uh, res)
    print("Mass of solution (average removed):", uh.integrate()[0])
    print("residual:", math.sqrt( res.as_petsc.dot(res.as_petsc) ) )

    error = math.sqrt(integrate(surface,(uh[0]-exact[0])**2,order=5)[0])

    l2error2 = inner(uh-exact, uh-exact)
    l2error2_gf = create.function("ufl", surface, "error", 5, l2error2)
    surface.writeVTK("laplace-surf", pointdata=[uh, l2error2_gf])
    return error

for i in range(4):
    error = compute()
    if i > 0:
        eoc = math.log(error/old_error) / math.log(0.5)
        print("size:", surface.size(0), "L2-error/eoc:",error,eoc)
    old_error = error
    surface.hierarchicalGrid.globalRefine(2)
