from __future__ import print_function

from mpi4py import MPI

import math
import ufl
import ufl.algorithms

import functools

import dune.models.elliptic
import dune.ufl
import dune.fem

import random

import abc

#base class for timestepping and setup of the phase field model

class phase_base(object):

    def __init__(self,sharpParam,auxParam):

        self.dimRange = 2
        self.dimDomain = 2

        #initializ grid variables
        self.initialGlobalRefine = auxParam[0]
        self.maxLevel = auxParam[1]
        self.saveInterval = auxParam[2]
        self.filename = auxParam[3]
        self.gridname = auxParam[4]

        #setup latent heat and anisotropic term
        self.LHeat = sharpParam[0]
        self.aniso = sharpParam[1]


        #set the initial count and time
        self.count = 0.
        self.t = 0.0

    #setup the phase field model from the sharp interface version
    def setupPhase(self):

        fac =  self.Gamma(self.Theta(self.un))
        tau = self.Tau(self.Theta(self.un))*self.Gamma(self.Theta(self.un))*self.eps_const*self.eps_const

        facdash = self.GammaDiff(self.Theta(self.un))

        #this is the interpolation function
        #can also come in front of K in temperature equations however here using phi_t
        pdash = 6*self.u[0]*(1-self.u[0])

        # right hand sie, time derivative part + explicit forcing in v
        a_ex = (ufl.inner(self.un, self.v) - self.LHeat * ufl.inner(self.un[0], self.v[1])) * ufl.dx
        # left hand side, heat equation in first variable + backward Euler in time


        #set up the diffusion tensor
        diag       = fac * fac
        offdiag    = -fac * facdash
        d0         = ufl.as_vector([diag, offdiag])
        d1         = ufl.as_vector([-offdiag, diag])

        #including the L(u) term here as well

        #6sqrt2) is the normalisation condtion to make sure \int \phi_{0z}^2 = 6sqrt(2) can be multiplied out at the end of the gibbs thompson condition
        s   = ufl.as_vector([self.dt_const / tau * self.u[0] * (1.0 - self.u[0]) * (self.u[0]-0.5) - self.dt_const /tau * self.eps_const / (6*ufl.sqrt(2)) * self.L(self.u[1])* pdash , self.LHeat * self.u[0]])

        a_im = (self.aniso* self.eps_const * self.eps_const * self.dt_const / tau * (ufl.inner(ufl.dot(d0, ufl.grad(self.u[0])), ufl.grad(self.v[0])[0]) + ufl.inner(ufl.dot(d1, ufl.grad(self.u[0])), ufl.grad(self.v[0])[1])) - self.dt_const * self.Sharp(self.u[1],self.v[1],self.time_const) + ufl.inner(self.u,self.v) - ufl.inner(s,self.v))  * ufl.dx

        #the end bit here is the sharp interface
        return [a_im, a_ex]


    #write current state to a file
    def write(self):
        from dune.models.elliptic import SourceWriter
        writer = SourceWriter("mymodel.hh")
        writer.openNameSpace('demo')
        model = dune.models.elliptic.compileUFL(a_im == a_ex, tempVars = False)
        model.write(writer, "MyModel")
        writer.closeNameSpace('demo')
        exit(1)


    def mark(self,element):
        solutionLocal = self.solution.localFunction(element)
        grad = solutionLocal.jacobian(element.geometry.domain.center)
        if grad[0].infinity_norm > 1.0:
          return self.hgrid.marker.refine if element.level < self.maxLevel else self.hgrid.marker.keep
        else:
          return self.hgrid.marker.coarsen


    # initial grid refinement
    def initialRefine(self):
        self.hgrid = self.grid.hierarchicalGrid

        #number of initial refinements to carry out on the grid
        self.grid.globalRefine(self.initialGlobalRefine)
        for i in range(0,self.maxLevel):
            self.hgrid.mark(lambda x: self.mark(x))
            self.hgrid.adapt([self.solution])
            self.hgrid.loadBalance([self.solution])
            self.solution.interpolate(self.initial_gf)
        return

    #setup all the spaces etc for the refinement
    def resetGrid(self):
        #setup function space etc
        self.uflSpace       = dune.ufl.Space(self.dimDomain, self.dimRange)
        self.uflScalarSpace = dune.ufl.Space(self.dimDomain, 1)
        self.u = ufl.TrialFunction(self.uflSpace)
        self.v = ufl.TestFunction(self.uflSpace)
        self.un = ufl.Coefficient(self.uflSpace)
        self.noise = ufl.Coefficient(self.uflSpace)

        #setup the time constant and dt so the ufl can be changed with each time step
        self.time_const = ufl.Constant(ufl.triangle)
        self.dt_const = ufl.Constant(ufl.triangle)
        self.eps_const = ufl.Constant(ufl.triangle)

        # basic setup
        # -----------
        self.grid       = dune.fem.leafGrid(self.gridname, "ALUSimplexGrid", dimgrid=self.dimDomain, refinement="conforming")
        self.spc        = dune.fem.create.space("Lagrange", self.grid, dimrange=self.dimRange, polorder=1)
        self.initial_gf = self.grid.globalGridFunction("initial", self.initial)
        self.noise_gf   = self.grid.globalGridFunction("noise", self.globalNoise)
        self.solution   = self.spc.interpolate(self.initial_gf, name="solution")
        self.solution_n = self.spc.interpolate(self.initial_gf, name="solution_n")
        self.noise_h    = self.spc.interpolate(self.noise_gf, name="noise_n")

        [a_im, a_ex] = self.setupPhase()

        # setup scheme
        self.model  = dune.models.elliptic.importModel(self.grid, dune.models.elliptic.compileUFL(a_im == a_ex)).get()
        self.scheme = dune.fem.create.scheme("FemScheme", self.solution, self.model, "scheme")
        self.model.setCoefficient(self.un, self.solution_n)

        # time loop setup
        self.count    = 0
        self.t        = 0.0
        self.saveStep = self.saveInterval
        self.vtk = self.grid.writeVTK(self.filename, pointdata=[self.solution], celldata=[self.grid.levelFunction(), self.grid.partitionFunction()], number=self.count)

        #carry out initial refinement on grid
        self.initialRefine()

        #initialize the time constant in ufl form

        try:
            self.model.setConstant(self.time_const, [0.])
        except:
            pass

        return 0

    def l2error(self,en,x):
        #what was solution here was initially uh
        y = en.geometry.position(x)
        exact_t =  lambda x: self.exact(self.t,x)
        val = self.solution.localFunction(en).evaluate(x) - exact_t(y)
        return [ math.sqrt( val[1]*val[1]) ];


    def getError(self):
        l2error_gf =  self.grid.localGridFunction( "error", self.l2error )
        error = self.grid.l2Norm(l2error_gf)
        return error

    def globalNoise(self,x):
        return [ 0, 0]

    def seteps(self,new_eps):
        self.eps = new_eps
        self.model.setConstant(self.eps_const,[new_eps])

    def setdt(self,new_dt):
        self.dt = new_dt
        self.model.setConstant(self.dt_const,[new_dt])

    #method to timestep once forward
    def nextTime(self):
        self.noise_h.interpolate(self.noise_gf)
        self.solution_n.assign(self.solution)
        self.scheme.solve(target=self.solution)
        self.t += self.dt

        #should add a try in here because if it isn't included it will produce error
        #change the time constant in the ufl form
        try:
            self.model.setConstant(self.time_const,[self.t])
        except:
            pass

        print('count: ',self.count,"t = ",self.t)
        if self.t > self.saveStep:
            self.saveStep += self.saveInterval
            self.count += 1
            self.vtk.write(self.filename, self.count)
        self.hgrid.mark(lambda x: self.mark(x))
        self.hgrid.adapt([self.solution])
        self.hgrid.loadBalance([self.solution])
        return 0


    #these are the specialised function that should be in the class one level up overwriteable

    #make Gamma Gamma_diff and initial abstract methods
    #here x is actually theta

    def Tau(self,theta):
        return 1.0/self.Gamma(theta)

    #L(u) = u in most situations
    def L(self,u):
        return u

    #here x is the solution function and y is the test function
    def Sharp(x,y,t):
        return ufl.inner(ufl.grad(x), ufl.grad(y))

    def Theta(self,phi):
        return ufl.atan_2(ufl.grad(phi[0])[1],ufl.grad(phi[0])[0])

    #exact solution of u
    def exact_u(self,t,x):
        return 0

    #exact solution of u and the phasefield
    def exact(self,t,x):
        #x[1] is the loc
        return [0.,exact_u(self,t,x)  ]


    #initial location of the interface and temperature 1 inside the region and 0 outside
    @abc.abstractmethod
    def initial(self,x):
        return

    @abc.abstractmethod
    def Gamma(self,theta):
        '''
        The anisotropic term written in terms of theta
        '''
        return

    @abc.abstractmethod
    def GammaDiff(self,theta):
        '''
        The differential of the gamma function defined above must be done by hand at the moment
        '''
        return
