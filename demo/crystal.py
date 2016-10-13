from __future__ import print_function

from mpi4py import MPI

import math
import ufl

import dune.ufl
import dune.common as common
import dune.fem as fem
import dune.create as create
from dune.fem.function import levelFunction, partitionFunction

def plot(solution,block=True):
    try:
        from matplotlib import pyplot
        from numpy import amin, amax, linspace
        pyplot.ion()

        # wouldn't it be better to have a  method
        # triagulation, data = solution.plot(0)
        # otherwise there is a constistency problem with the level parameter
        triangulation = solution.grid.triangulation(0)
        data = solution.pointData(0)

        levels = linspace(amin(data[:,0]), amax(data[:,0]), 256)

        pyplot.gca().set_aspect('equal')
        pyplot.triplot(solution.grid.triangulation(), antialiased=True, linewidth=0.2, color='black')
        pyplot.tricontourf(triangulation, data[:,0], cmap=pyplot.cm.rainbow, levels=levels)
        pyplot.pause(0.005)
        if block:
            pyplot.show(block=True)
    except ImportError:
        pass

def compute():
    # problem parameters
    # ------------------
    dimRange     = 2
    dimDomain    = 2
    maxLevel     = 10
    dt           = 5.e-4
    endTime      = 0.2
    saveinterval = 0.01
    order        = 1

    ## model taken from www.ctcms.nist.gov/fipy/examples/phase/generated/examples.phase.anisotropy.html
    alpha        = 0.015
    tau          = 3.e-4
    kappa1       = 0.9
    kappa2       = 20.
    c            = 0.02
    N            = 6.
    def initial(x):
        r  = (x-[6,6]).two_norm
        return [ 0 if r>0.3 else 1, -0.5 ]

    # set up left and right hand side models
    # --------------------------------------
    uflSpace = dune.ufl.Space(dimDomain, dimRange)
    u = ufl.TrialFunction(uflSpace)
    v = ufl.TestFunction(uflSpace)
    un = ufl.Coefficient(uflSpace)

    # right hand sie (time derivative part + explicit forcing in v)
    a_ex = (ufl.inner(un, v) - ufl.inner(un[0], v[1])) * ufl.dx
    # left hand side (heat equation in first variable + backward Euler in time)
    psi        = ufl.pi/8.0 + ufl.atan_2(ufl.grad(un[0])[1], (ufl.grad(un[0])[0]))
    Phi        = ufl.tan(N / 2.0 * psi)
    beta       = (1.0 - Phi*Phi) / (1.0 + Phi*Phi)
    dbeta_dPhi = -2.0 * N * Phi / (1.0 + Phi*Phi)
    fac        = 1.0 + c * beta
    diag       = fac * fac
    offdiag    = -fac * c * dbeta_dPhi
    d0         = ufl.as_vector([diag, offdiag])
    d1         = ufl.as_vector([-offdiag, diag])
    m          = u[0] - 0.5 - kappa1 / ufl.pi*ufl.atan(kappa2*u[1])
    s          = ufl.as_vector([dt / tau * u[0] * (1.0 - u[0]) * m, u[0]])

    a_im = (alpha*alpha*dt / tau * (ufl.inner(ufl.dot(d0, ufl.grad(u[0])), ufl.grad(v[0])[0]) + ufl.inner(ufl.dot(d1, ufl.grad(u[0])), ufl.grad(v[0])[1]))
           + 2.25 * dt * ufl.inner(ufl.grad(u[1]), ufl.grad(v[1])) + ufl.inner(u,v) - ufl.inner(s,v)) * ufl.dx

    # basic setup
    # -----------
    grid       = create.view("adaptive", create.grid("ALUConform", "../data/crystal-2d.dgf", dimgrid=dimDomain))
    spc        = create.space("Lagrange", grid, dimrange=dimRange, order=order)
    initial_gf = create.function("global", grid, "initial", order+1, initial)
    solution   = spc.interpolate(initial_gf, name="solution")
    solution_n = spc.interpolate(initial_gf, name="solution_n")

    # setup scheme
    # ------------
    model  = create.model("elliptic", grid, a_im == a_ex, coefficients={un:solution_n} )
    scheme = create.scheme("h1", solution, model,
            parameters={
            "fem.solver.newton.tolerance": 1e-5,
            "fem.solver.newton.linabstol": 1e-8,
            "fem.solver.newton.linreduction": 1e-8,
            "fem.solver.newton.verbose": 1,
            "fem.solver.newton.linear.verbose": 1}\
            )


    # marking strategy
    # ----------------
    def mark(element):
        marker = common.Marker
        solutionLocal = solution.localFunction(element)
        grad = solutionLocal.jacobian(element.geometry.domain.center)
        if grad[0].infinity_norm > 1.0:
          return marker.refine if element.level < maxLevel else marker.keep
        else:
          return marker.coarsen

    # initial grid refinement
    # -----------------------
    hgrid = grid.hierarchicalGrid
    hgrid.globalRefine(2)
    for i in range(0,maxLevel):
        hgrid.mark(mark)
        fem.adapt(hgrid, [solution])
        fem.loadBalance(hgrid, [solution])
        solution.interpolate(initial_gf)

    # time loop
    # ---------
    count    = 0
    t        = 0.0
    savestep = saveinterval
    vtk = grid.writeVTK("crystal",
            pointdata=[solution],
            celldata=[levelFunction(grid), partitionFunction(grid)],
            number=count)

    while t < endTime:
        solution_n.assign(solution)
        scheme.solve(target=solution)
        t += dt
        print('count: ',count,"t = ",t)
        if t > savestep:
            savestep += saveinterval
            count += 1
            vtk.write("crystal", count)
            plot(solution, False)
        hgrid.mark(mark)
        fem.adapt(hgrid, [solution])
        fem.loadBalance(hgrid, [solution])

    plot(solution, True)
    print("END")

compute()
