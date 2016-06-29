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
x = model.spatialCoordinate()
a = (inner(u,v) + inner(grad(u),grad(v)))*dx(0)
c0 = cos(2*math.pi*x[0])
c1 = cos(2*math.pi*x[1])
g = [c0*c1,c0*c0]

# model.generate(a,exact=exact)
model.generate(a,diric={1:g},exact=exact)
model.write(exact=exact, name="LaplaceMass")
