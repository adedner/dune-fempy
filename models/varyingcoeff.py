print 'Warning: this model includes a coefficient for the mass and diffusion coefficients -\
       this can not yet be set so the model does not compile'
from dune.models.femufl import *

model    = DuneUFLModel(2,1,'VaryingCoeff') # this is utility.init and sets the dim range

#########################################
# define exact solution for testing
# as sympy expression
#########################################
exact = [sympy.cos(2*math.pi*model.x0)*sympy.cos(2*math.pi*model.x1)]
#########################################
# define main ufl form a
#########################################
u = model.trialFunction()
v = model.testFunction()
x = model.spatialCoordinate()
D = model.coefficient('diffusion')
m = model.coefficient('mass')
dx0 = dx(0)
a = ( m[0]*inner(u,v) + D[0]*inner(grad(u),grad(v)) )*dx0
L = sin(x[0])*sin(x[1])*v[0]*dx(0)
model.generate(a,L,exact)
