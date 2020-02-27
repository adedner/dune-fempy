import math, os, sys, time

import dune.create as create
from dune.grid import cartesianDomain
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction, replace,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
import dune.fem
from math import pi,log10

from dune.fem.function import integrate

from ufl import cos, sin, exp, sqrt

from dune.fem import parameter

from .presscorrectclass import PressCorrect
from .pressupdateclass import PressUpdate
from .computevelclass import ComputeVel
from .combspcmtds.altpresscorrectclass import AltPressCorrect
from .combspcmtds.altcomputevelclass import AltComputeVel

class Scheme:
    chorin = 0
    ipcs = 1
    ripcs = 2
    hdivripcs =3

    names = {chorin:"chorin",
             ipcs:"ipcs",
             ripcs:"ripcs",
             hdivripcs:"hdivripcs"}

    @staticmethod
    def toString(method):
        return Scheme.names[method]
    def fromString(name):
        return next( (method for method, n in Scheme.names.items() if n == name), None)

    def __init__(self, method, model, solver=None, parameters=None):
        self.method = method
        if model.useCombSPACE == True:
            self.computeTentativeVel  = AltComputeVel(model,method,velCorrectionStep=False)
            self.pressureCorrection   = AltPressCorrect(model,method)
            self.velupdt              = AltComputeVel(model,method,velCorrectionStep=True)
        else:
            self.computeTentativeVel  = ComputeVel(model,method,velCorrectionStep=False)
            self.pressureCorrection   = PressCorrect(model,method)
            self.velupdt              = ComputeVel(model,method,velCorrectionStep=True)

        #
        self.theta = 1
        self.alpha = 0.5
        self.beta  = 1-self.alpha
        self.deltaT = model.deltaT
    def __str__(self):
        return Scheme.toString(self.method)
    def solve(self,simTime, target):
        print( 'Solve step 1 - ComputeVel' )
        if self.method  == ComputeVel.hdivripcs:
            self.computeTentativeVel.l2project(target)
        self.computeTentativeVel.prepare(simTime+self.deltaT,1./(self.theta*self.deltaT),self.alpha,target,target)
        self.velupdt.saveold(simTime+self.deltaT,1/(self.deltaT),self.alpha,target,target)
        self.computeTentativeVel.solve(target)
        print( 'Solve step 2 - Pressure correction' )
        self.pressureCorrection.saveold(simTime+self.deltaT,1/(self.deltaT),self.alpha,target,target)
        self.pressureCorrection.prepare(simTime+self.deltaT,1/(self.deltaT),self.alpha,target,target)
        self.pressureCorrection.solve(target)
        if self.method  == ComputeVel.hdivripcs:
            print( 'Solve step 3 - Hdiv projection' )
            self.velupdt.hdivproject(target)
        else:
            print( 'Solve step 3 - Velocity correction' )
            self.velupdt.prepare(simTime+self.deltaT,1./(self.deltaT),self.alpha,target,target)
            self.velupdt.solve(target)