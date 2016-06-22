print 'Warning: this model does not lead a converging scheme'
from dune.models.femufl import *

model    = DuneUFLModel(2,3,'System') # this is utility.init and sets the dim range

#########################################
# define exact solution for testing
# as sympy expression
#########################################
exact = [sympy.cos(2*math.pi*model.x0)*sympy.cos(2*math.pi*model.x1),
         model.x0*model.x1,
         (1+model.x0*model.x0)
        ]
#########################################
# define main ufl form a
#########################################
u = model.trialFunction()
v = model.testFunction()
x = model.spatialCoordinate()
a = ( inner(u,u)*u[0]/3*v[0] + (inner(u,u)+2)*dot(grad(u[0]), grad(v[0])) ) * dx(0)
a = a + dot(grad(u[1]),grad(v[1])) * dx(0)
a = a + u[2]*v[2] * dx(0)

model.generate(a,exact=exact)
