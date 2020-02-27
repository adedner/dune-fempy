from dune import create
from dune.fem.function import integrate

from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
from ufl import cos, sin, exp, sqrt
from math import pi,log10
from ufl import *

class PressUpdate:
    chorin = 0
    ipcs = 1
    css = 2
    names = {chorin:"chorin",
             ipcs:"ipcs",
             css:"css"}
    def __init__(self,model,method):
        spc  = model.spc
        Re       = model.Re
        forcing  = model.f
        #velocity bnd conditions
        bcs      = model.bcs_press
        spcU = spc[0]
        spcP = spc[1]
        u  = TrialFunction(spcU)
        v  = TestFunction(spcU)
        p  = TrialFunction(spcP)
        q  = TestFunction(spcP)
        mu = NamedConstant(spcU, "mu")
        nu = NamedConstant(spcU, "nu")
        x     = SpatialCoordinate(spcU)
        time  = NamedConstant(spcU, "time")     # current time
        self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")
        self.altoldVelo = spcU.interpolate([0,]*spcU.dimRange,name="altoldU")

        self.oldPress = spcP.interpolate([0],name="oldP")
        self.altoldPress = spcP.interpolate([0],name="altoldP")


        Du       = 0.5 *(grad(u)+grad(u).T)
        oldDu    = 0.5 *(grad(self.oldVelo)+grad(self.oldVelo).T)
        norm2_Du = (inner(Du, Du))

        nu    = 1./Re
        gamma_t = lambda t: exp(-2.*pi*pi*nu*t)
        exact_u = as_vector( [ -1*cos(1.*pi*x[0])*sin(1.*pi*x[1])*gamma_t(time), sin(1.*pi*x[0])*cos(1.*pi*x[1])*exp(-2.*pi*pi*nu*time) ] )
        exact_p = as_vector( [ -0.25*(cos(2.*pi*x[0])+cos(2.*pi*x[1]))*exp(-4.*pi*pi*nu*time) ] )

        #if method == ComputeVel.css:
        if model.order == 1:
            ps = self.oldPress[0]
        else:
            ps = 2*self.oldPress[0]- self.altoldPress[0]
        f = inner(ps,q[0]) *dx  + inner(self.oldPress[0],q[0]) *dx - ((1/Re)*q[0]*div(self.oldVelo) ) * dx
        #LHS
        pressCorrectmodel = inner(p[0],q[0])  * dx

        solverParameter={"fem.solver.newton.absolutetol": 1e-11,
                         "fem.solver.newton.reductiontol": 1e-11,
                         "fem.solver.newton.tolerance": 1e-10,
                         "fem.solver.newton.verbose": "true",
                         "fem.solver.newton.linear.verbose": "false"}
        # bcs_press = [DirichletBC(spcP,[0],1)]

        #use direct implementation of exact sol. only in case of time dependant boundary conditions
        self.scheme = create.scheme("h1", [pressCorrectmodel==f,DirichletBC(spcP,[exact_p[0]],1)], spcP, solver="gmres", parameters = solverParameter)
        # self.scheme = create.scheme("h1", [pressCorrectmodel==f,*bcs], spcP, solver="gmres", parameters = solverParameter)


    def prepare(self,t,mu,nu,oldSolution,oldSolution1):
        try:
            self.scheme.model.time = t
        except: pass
        #self.scheme.model.mu = mu
        # self.scheme.model.nu = nu
        self.oldVelo.assign(oldSolution[0])
        self.oldVelo.assign(oldSolution1[0])
        self.oldPress.assign(oldSolution[1])
        self.altoldPress.assign(oldSolution1[1])


    def solve(self,target):
        self.scheme.solve(target=target[1])
