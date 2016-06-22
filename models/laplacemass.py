from dune.models.femufl import *

model    = DuneUFLModel(2,2,'LaplaceMass') # this is utility.init and sets the dim range

#########################################
# define exact solution for testing
# as sympy expression
#########################################
c0 = sympy.cos(2*math.pi*model.x0)
c1 = sympy.cos(2*math.pi*model.x1)
exact = [c0*c1,c0*c1]
# exact = [sympy.cos(2*math.pi*model.x0)*sympy.cos(2*math.pi*model.x1)]
#########################################
# define main ufl form a
#########################################
u = model.trialFunction()
v = model.testFunction()
x = model.spatialCoordinate()
dx0 = dx(0)
a = (inner(u,v) + inner(grad(u),grad(v)))*dx0

model.generate(a,exact=exact)
