from dune import create
from dune.fem.function import integrate


class Burgers:
    def __init__(self,spc,forcing):
        spcU = spc[0]
        spcP = spc[1]
        mu = NamedConstant(spcU, "mu")
        nu = NamedConstant(spcU, "nu")
        self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")
        self.oldPress = spcU.interpolate([0]*spcP.dimRange,name="oldP")

        f = (dot(forcing,v) +
             mu*dot(self.oldVelo,v) - theta2 * deltaT * alpha * nu*inner(grad(self.oldVelo)+grad(self.oldVelo).T, grad(v)) -
             inner(grad(self.oldVelo), outer(v, self.oldVelo)) -1* theta2 *deltaT * inner(self.oldPress,div(v)) ) * dx

        old_solution = velocity.copy();
        old_pressure = pressure.copy();

        old_solution.name = "uOld"
        u_n = old_solution
        old_pressure.name = "pOld"
        p_n = old_pressure
        # a = (inner(u - u_n, v) + inner(grad(u), outer(v, u_n))) * dx
        a = (inner(u - u_n, v) + theta2 * deltaT * beta * nu * inner(grad(u)+grad(u).T, grad(v)) + theta2 * deltaT * inner(grad(u), outer(v, u))) * dx
        a += ( theta2 * deltaT * alpha * nu * inner(grad(u_n)+grad(u_n).T, grad(v))) * dx
        a += (-1* theta2 *deltaT * inner(p_n[0],div(v)) )* dx


        model = create.model("integrands", grid, a == f)

        solverParameter={"fem.solver.newton.linabstol": 1e-11,
                         "fem.solver.newton.linreduction": 1e-11,
                         "fem.solver.newton.tolerance": 1e-10,
                         "fem.solver.newton.verbose": "true",
                         "fem.solver.newton.linear.verbose": "false"}

        self.scheme = create.scheme("galerkin", a == 0, spcU, solver="gmres", parameters = solverParameter)


    def prepare(self,t,mu,nu,velo):
        try:
            self.scheme.time = t
        except: pass
            old_solution.assign(velo)
            old_pressure.assign(velo)

    def solve(self,target):
        velocity = target[0]
        pressure = target[1]
        scheme.solve(target=velocity)
