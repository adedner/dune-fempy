from dune.models.femufl import *

model    = DuneUFLModel(2,1,'NonLinear') # this is utility.init and sets the dim range

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
a = ( u[0]*u[0]*u[0]/3*v[0] + (u[0]*u[0]+2)*dot(grad(u[0]), grad(v[0])) ) * dx(0)

model.generate(a,exact=exact)
