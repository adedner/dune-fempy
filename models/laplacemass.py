from dune.models.femufl import *

model    = DuneUFLModel(2,2) # this is utility.init and sets the dim range

#########################################
# define exact solution for testing
# as sympy expression
#########################################
c0 = sympy.cos(2*math.pi*model.x0)
c1 = sympy.cos(2*math.pi*model.x1)
exact = [c0*c1,c0*c0]
#########################################
# define main ufl form a
#########################################
u = model.trialFunction()
v = model.testFunction()
a = (inner(u,v) + inner(grad(u),grad(v)))*dx(0)

model.generate(a, exact=exact)
model.write(exact=exact, name="LaplaceMass")
