from dune import create
from dune.fem.function import integrate

from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant

class Burgers:
    def __init__(self,spc,forcing,bcs,Re):
        spcU = spc[0]
        spcP = spc[1]
        u  = TrialFunction(spcU)
        v  = TestFunction(spcU)
        p  = TrialFunction(spcP)
        q  = TestFunction(spcP)
        mu = NamedConstant(spcU, "mu")
        nu = NamedConstant(spcU, "nu")
        self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")
        self.oldPress = spcP.interpolate([0],name="oldP")

        f = (dot(forcing,v) +
             mu*dot(self.oldVelo,v) - nu/Re*inner(grad(self.oldVelo)+grad(self.oldVelo).T, grad(v)) -
             self.oldPress[0]*div(v) ) * dx

        a = (mu*inner(u, v) + (1-nu)/Re * inner(grad(u)+grad(u).T, grad(v)) + inner(grad(u), outer(v, u))) * dx

        solverParameter={"fem.solver.newton.absolutetol": 1e-11,
                         "fem.solver.newton.reductiontol": 1e-11,
                         "fem.solver.newton.tolerance": 1e-10,
                         "fem.solver.newton.verbose": "true",
                         "fem.solver.newton.linear.verbose": "false"}
        self.scheme = create.scheme("h1", [a==f,*bcs], spcU, solver="gmres", parameters = solverParameter)


    def prepare(self,t,mu,nu,oldSolution):
        try:
            self.scheme.time = t
        except: pass
        self.scheme.model.mu = mu
        self.scheme.model.nu = nu
        self.oldVelo.assign(oldSolution[0])
        self.oldPress.assign(oldSolution[1])

    def solve(self,target):
        self.scheme.solve(target=target[0])
