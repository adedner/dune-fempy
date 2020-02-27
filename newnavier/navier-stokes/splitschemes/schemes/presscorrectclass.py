from dune import create
from dune.fem.function import integrate

from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
from ufl import cos, sin, exp, sqrt
from math import pi,log10
from ufl import *

class PressCorrect:
    chorin = 0
    ipcs = 1
    ripcs = 2
    hdivripcs =3

    names = {chorin:"chorin",
             ipcs:"ipcs",
             ripcs:"ripcs",
             hdivripcs:"hdivripcs"}

    def __init__(self,model,method):
        spc  = model.spc
        Re       = model.Re
        forcing  = model.f
        bcs_press      = model.bcs_press
        self.method = method

        spcU = spc[0]
        spcP = spc[1]
        u  = TrialFunction(spcU)
        v  = TestFunction(spcU)
        p  = TrialFunction(spcP)
        q  = TestFunction(spcP)
        mu = NamedConstant(spcU, "mu")
        nu = NamedConstant(spcU, "nu")
        x     = SpatialCoordinate(spcU)
        time  = NamedConstant(spcU, "t")
        n = FacetNormal(spcU.cell())


        self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")
        self.prevoldVelo = spcU.interpolate([0,]*spcU.dimRange,name="altoldU")

        self.oldPress = spcP.interpolate([0],name="oldPress")
        self.prevoldPress = spcP.interpolate([0],name="prevoldPress")

        Du       = 0.5 *(grad(u)+grad(u).T)
        oldDu    = 0.5 *(grad(self.oldVelo)+grad(self.oldVelo).T)
        norm2_Du = (inner(Du, Du))
        nu    = 1./Re

        if self.method == PressCorrect.ipcs or self.method == PressCorrect.ripcs or self.method == PressCorrect.hdivripcs:
            f = (-1*mu*q[0]*div(self.oldVelo) ) * dx
            pressCorrectmodel =  inner(grad(p[0]-self.oldPress[0]), grad(q[0])) * dx

            if self.method == PressCorrect.ripcs:
                pressCorrectmodel +=  inner(grad(nu*div(self.oldVelo)), grad(q[0])) * dx

            # if model.useDGPress == True:
            #     penal = 7.5 * 1006
            #
            #     # a =  1/Re*2*inner(grad(u), grad(v)) * dx
            #     a =  1/Re*2*(inner(outer(jump(p[0]), n('+')), avg(grad(q[0]))) + inner(avg(grad(p[0])), outer(jump(q[0]), n('+')))) * dS
            #     a += 1/Re*2* penal * inner(jump(p[0]), jump(q[0])) * dS
            #     a -=  1/Re*2*(inner(outer(p[0]-exact_p[0], n), grad(q[0])) + inner(grad(p[0]-exact_p[0]), outer(q[0], n))) * ds
            #     a +=  1/Re*2* penal * inner(p[0]-exact_p[0], q[0]) * ds
            #     pressCorrectmodel  +=  a
        elif self.method == PressCorrect.chorin:
            f = (-1*mu*q[0]*div(self.oldVelo) ) * dx

            pressCorrectmodel =  inner(grad(p[0]), grad(q[0])) * dx
            # if model.useDGPress == True:
            #     penal = 7.5 * 1006
            #
            #     # a =  1/Re*2*inner(grad(u), grad(v)) * dx
            #     a =  1/Re*2*(inner(outer(jump(p[0]), n('+')), avg(grad(q[0]))) + inner(avg(grad(p[0])), outer(jump(q[0]), n('+')))) * dS
            #     a += 1/Re*2* penal * inner(jump(p[0]), jump(q[0])) * dS
            #     a -=  1/Re*2*(inner(outer(p[0]-exact_p[0], n), grad(q[0])) + inner(grad(p[0]-exact_p[0]), outer(q[0], n))) * ds
            #     a +=  1/Re*2* penal * inner(p[0]-exact_p[0], q[0]) * ds
            #     pressCorrectmodel  +=  a

        if model.useDGPress == True:
            self.scheme = create.scheme("galerkin", [pressCorrectmodel==f], spcP, solver="cg")
        else:
            self.scheme = create.scheme("galerkin", [pressCorrectmodel==f,*bcs_press], spcP, solver="cg")

    def prepare(self,t,mu,nu,oldSolution,oldSolution1):
        self.scheme.model.t = t
        self.scheme.model.mu = mu
        self.oldVelo.assign(oldSolution[0])
        self.oldPress.assign(oldSolution[1])

    def saveold(self,t,mu,nu,oldSolution,oldSolution1):
        self.prevoldPress.assign(oldSolution[1])

    def solve(self,target):
        self.scheme.solve(target=target[1])
