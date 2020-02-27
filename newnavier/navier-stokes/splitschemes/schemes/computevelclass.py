from dune import create
from dune.fem.function import integrate
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx,ds, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
from ufl import cos, sin, exp, sqrt
from math import pi,log10
from ufl import *

class ComputeVel:
    chorin = 0
    ipcs = 1
    ripcs = 2
    hdivripcs =3

    names = {chorin:"chorin",
             ipcs:"ipcs",
             ripcs:"ripcs",
             hdivripcs:"hdivripcs"}

    def epsilon(u):
        return 0.5*(grad(u) + grad(u).T)
    def sigma(u, p, nu):
        return 2*nu*epsilon(u) - p*Identity(u.cell().d)

    def __init__(self,model,method,velCorrectionStep=False):
        spc      = model.spc
        Re       = model.Re
        forcing  = model.f
        bcs      = model.bcs
        self.method = method

        spcU = spc[0]
        spcP = spc[1]
        self.spcU = spcU
        dimension = spcU.grid.dimension
        u  = TrialFunction(spcU)
        v  = TestFunction(spcU)
        p  = TrialFunction(spcP)
        q  = TestFunction(spcP)
        mu = NamedConstant(spcU, "mu")
        x     = SpatialCoordinate(spcU)
        time  = NamedConstant(spcU, "t")
        n = FacetNormal(spcU.cell())

        self.rtspc = create.space("raviartThomas", model.grid, dimrange=2, order=2, storage="istl")

        self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")
        self.oldPress = spcP.interpolate([0],name="oldP")
        self.prevoldPress = spcP.interpolate([0],name="altoldP")

        Du       = 0.5 *(grad(u)+grad(u).T)
        Dv       = 0.5 *(grad(v)+grad(v).T)

        oldDu    = 0.5 *(grad(self.oldVelo)+grad(self.oldVelo).T)
        norm2_Du = (inner(Du, Du))
        nu    = 1./Re

        UU      =  0.5 * (u+self.oldVelo)
        p0      =  self.oldPress[0]
        # p0      =  self.oldPress[0] + self.prevoldPress[0]
        pbar    = 0.5*(p0 + p[0])
        Sig     =  2*(1/Re)*0.5 *(grad(UU)+grad(UU).T)

        if forcing is not None:
            f = inner(forcing,v)*dx
        else:
            f = 0


        if self.method == ComputeVel.ipcs or self.method == ComputeVel.ripcs or self.method == ComputeVel.hdivripcs:
            if velCorrectionStep:
                mainModel = inner(u,v)*dx
                rhs = inner(self.oldVelo,v)*dx  - ((1/mu)*inner(v, grad(self.oldPress[0]-self.prevoldPress[0]))) * dx

                if self.method == ComputeVel.ripcs:
                    rhs +=  - (1/mu)*inner(v, grad(nu*div(self.oldVelo))) * dx
                #LHS
            else:
                rhs  =  f + (mu*inner(self.oldVelo,v) - inner(grad(self.oldVelo), outer(v, self.oldVelo))) * dx
                rhs +=  inner(Dv, p0*Identity(dimension)) * dx
                rhs +=  -1 * inner(dot(v, n), p0)*ds
                #LHS
                #set to 0 for periodic bnd cond otherwise 1
                betaflag = 1
                mainModel  = (mu*inner(u,v) + inner(Dv, Sig)) * dx - betaflag*nu*inner(outer(v, n), grad(UU).T)*ds
                if model.useDG == True:
                    penal = 7.5 * 1006
                    a =  1/Re*2*(inner(outer(jump(UU), n('+')), avg(grad(v))) + inner(avg(grad(UU)), outer(jump(v), n('+')))) * dS
                    a += 1/Re*2* penal * inner(jump(UU), jump(v)) * dS
                    a -=  1/Re*2*(inner(outer(UU-exact_u, n), grad(v)) + inner(grad(UU-exact_u), outer(v, n))) * ds
                    a +=  1/Re*2* penal * inner(UU-exact_u, v) * ds
                    mainModel  +=  a
        elif self.method == ComputeVel.chorin:
            if velCorrectionStep:
                rhs = (inner(self.oldVelo,v) - (1/mu)*inner(v, grad(self.oldPress[0]))) * dx
                #LHS
                mainModel = inner(u,v)  * dx
            else:
                rhs= f + (mu*inner(self.oldVelo,v) - inner(grad(self.oldVelo), outer(v, self.oldVelo))) * dx
                #LHS
                abs_du = (inner(Du, Du))
                d = 0.00
                pnb = 2
                b = 1/Re*2*(pow(d**(1) + (abs_du+1.0e-10)**0.5, (pnb-2)/1)*inner(Du, grad(v)) ) * dx
                mainModel  = (mu*inner(u,v))*dx + b

                if model.useDG == True:
                    penal = 7.5 * 1006
                    a =  1/Re*2*(inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
                    a += 1/Re*2* penal * inner(jump(u), jump(v)) * dS
                    a -=  1/Re*2*(inner(outer(u-exact_u, n), grad(v)) + inner(grad(u-exact_u), outer(v, n))) * ds
                    a +=  1/Re*2* penal * inner(u-exact_u, v) * ds
                    mainModel  +=  a
        else:
            print("Invalid name for method:")

        solverParameter={}

        if model.useDG == True:
            self.scheme = create.scheme("galerkin", [mainModel==rhs], spcU, solver="cg", parameters = solverParameter)
        else:
            self.scheme = create.scheme("galerkin", [mainModel==rhs,*bcs], spcU, solver="cg", parameters = solverParameter)

    def prepare(self,time=0,mu=0,nu=1,oldSolution=None,oldSolution1=None):
        # check that variables exist
        # b = self.scheme.model.t

        self.scheme.model.t = time
        # print ("Scheme: time is set to ", self.scheme.model.time)
        self.scheme.model.mu = mu
        if oldSolution:
            self.oldVelo.assign(oldSolution[0])
            self.oldPress.assign(oldSolution[1])

    def saveold(self,t,mu,nu,oldSolution,oldSolution1):
        self.prevoldPress.assign(oldSolution[1])

    def solve(self,target):
        self.scheme.solve(target[0])

    def hdivproject(self,target):
        target[0] = self.rtspc.interpolate( target[0])
        # target[0].interpolate( target[0])

    def l2project(self,target):
        target[0] = self.spcU.interpolate( target[0])
