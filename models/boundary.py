from dune.models.femufl import *

model    = DuneUFLModel(2,2) # this is utility.init and sets the dim range

#########################################
# define main ufl form a
#########################################
u = model.trialFunction()
v = model.testFunction()
x = model.spatialCoordinate()
bnd = model.coefficient("bnd", 2)
g1 = [x[1]*cos(x[0]),sin(x[0])]
g2 = [bnd,0]
a = (inner(u,v) + inner(grad(u),grad(v)))*dx(0) + inner(u,u)*v[0]*ds(0)

model.generate(a, diric = {1:g1, 2:g2})
model.write(name="Boundary")
