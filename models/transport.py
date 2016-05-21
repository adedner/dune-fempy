from dune.models.femufl import *

model    = DuneUFLModel(2,1,'Transport') # this is utility.init and sets the dim range
vecmodel = DuneUFLModel(2,2,'Velocity') # this is utility.init and sets the dim range

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
velo = vecmodel.coefficient('velocity')
diff = model.coefficient('diffusion')
a = ( (dot(velo,grad(u[0]))+0.01*u[0])*v[0] + diff[0]*inner(grad(u[0]),grad(v[0])) ) * dx(0)
L = 10./(1.+(x[0]*x[0]+x[1]*x[1])**4 )  *  v[0]*dx(0)

velo = [-model.x1,model.x0]
# model.setCoefficient("velocity",velo)
model.setCoefficient("diffusion",[1+0.*model.x0*model.x0])

# model.generateFromExact(a,exact)
model.generate(a,L,exact)
