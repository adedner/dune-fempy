print 'Warning: this model does not lead a converging scheme - \
       possible issue with boundary conditions or forcing'
from dune.models.femufl import *

model    = DuneUFLModel(2,1) # this is utility.init and sets the dim range

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
a = ( inner(grad(u),grad(v))/sqrt(1+inner(grad(u),grad(u))) ) * dx(0)

model.generate(a,exact=exact)
model.write(exact=exact, name="MeanCurvature")
