from phase_base import *

#define class with all the user functions for specifically sharp interface model
class sharpModel(phase_base):
    def __init__(self,auxParam):

        LHeat = 1.      #latent heat in the jump condition at the interface
        aniso = 1.      #stregnth of anisotropic mean curvature

        #parameters directly relating the the sharp interface model
        sharpParam = (LHeat,aniso)

        #make sure initialize all variables in base class
        phase_base.__init__(self,sharpParam,auxParam)

        self.N = 6.

    def initial(self,x):
        r = (x-[6,6]).two_norm
        return [ 0 if r > 0.3 else 1, -0.5]

    def Theta(self,phi):
        return  ufl.tan(self.N / 2.0 * ( ufl.pi/8.0 +  ufl.atan_2(ufl.grad(phi[0])[1], ufl.grad(phi[0])[0])))

    def Tau(self,theta):
        alpha = 4./3
        return (1.0/self.Gamma(theta))*alpha

    def GammaDiff(self,theta):
        dbeta_dPhi =  -2.0 * self.N * theta / (1.0 + theta*theta)
        return 0.02*dbeta_dPhi

    @staticmethod
    def Gamma(theta):
        return 1.0+0.02* ( (1.0 - theta*theta) / (1.0 + theta*theta))

    @staticmethod
    def L(u):
        gamma = 27.00948948
        return gamma*ufl.atan(20.*u)

    #notice this also has a t dependency so that time dependent source functions can be included
    @staticmethod
    def Sharp(u,v,t):
        return -2.25*ufl.inner(ufl.grad(u), ufl.grad(v))

    #exact solution for calculating the error
    def exact_u(self,t,x):
        #  = (x-[6,6]).two_norm
        return 0

#computation parameters not to do with the sharp interface model defined outside the class
maxLevel = 12
initialGlobalRefine = 2                 #initial number of initial global refinements
saveInterval = 0.001                    #interval to save the files
filename = 'testingfile'                #name of file to saved
gridname = '../data/crystal-2d.dgf'     #location of initial grid file

#tuple of extra parameters for passing
auxParam = (initialGlobalRefine,maxLevel,saveInterval,filename,gridname)

#set up
CrystalModel = sharpModel(auxParam)

#set up function spaces and time loop resets all grids and time
CrystalModel.resetGrid()

#dt is set like this so that it can be changed at each time step if needed as is common in these types of problems. These are defined as constant types in ufl so eps can be changed without recompliing
CrystalModel.setdt(5.e-4)
CrystalModel.seteps(0.015)

#loop through each time
while CrystalModel.t < 0.2:
    CrystalModel.nextTime()

    #if defined an exact solution can get with model1.getError()
    #returns the l2 error
    #print('error=',  model1.getError())

print("END")
