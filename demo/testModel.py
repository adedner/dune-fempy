"""Compiles and solves different models (for code testing).
"""
############################################################################
# initialize
############################################################################
from __future__ import print_function
import sys,os,math, getopt
sys.path.append("../python")
import timeit

import sympy
import ufl

import dune.models.femufl as duneuflmodel
import dune.fem.grid as grid
import dune.fem.scheme as scheme

grid2d = grid.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)
grid2d.globalRefine(3)

############################################################################
# build a transport model -eps Laplace(u) + velocity.grad(u) + gamma u = f
# with unknown velocity field
############################################################################
model    = duneuflmodel.DuneUFLModel(2,1,'Transport')
vecmodel = duneuflmodel.DuneUFLModel(2,2,'Velocity')
exact = [sympy.cos(2*math.pi*model.x0)*sympy.cos(2*math.pi*model.x1)]
u = model.trialFunction()
v = model.testFunction()
x = model.spatialCoordinate()
velo = vecmodel.coefficient('velocity')
diff = model.coefficient('diffusion')
a = ( (ufl.dot(velo,ufl.grad(u[0]))+0.01*u[0])*v[0] + diff[0]*ufl.inner(ufl.grad(u[0]),ufl.grad(v[0])) ) * ufl.dx(0)
L = 10./(1.+(x[0]*x[0]+x[1]*x[1])**4 )  *  v[0]*ufl.dx(0)
model.setCoefficient("diffusion",[1+0.1*model.x0*model.x0])

start_time = timeit.default_timer()
model.generate(a,L,exact)
Model = model.makeAndImport(grid2d)
print("Building TransportModel took ", timeit.default_timer() - start_time, "s")
m = Model.get()

s = scheme.scheme("FemScheme", grid2d, m, "transport", polorder=1)

############################################################################
# Now use different way to solve this problem with given velocity field
############################################################################

print("define velocity field in global coordinates (using different approaches)")
expr1 = gf.MathExpression(["-(x1-0.5)","x0-1./2."])
expr2 = gf.SympyExpression(["-(x1-0.5)","x0-1./2."])
def expr_func(x,r):
    r[0] = -(x[1]-0.5)
    r[1] = (x[0]-0.5)
expr3 = gf.FuncExpression(2,expr_func)
velocityGlobal = grid2d.getGlobal("global_velocity",expr1)

m.setvelocity(velocityGlobal)
s.solve()
out_gl = grid2d.vtkOutput()
out_gl.add(s.solution())
print('solution: ', s.solution())
out_gl.add(velocityGlobal)
out_gl.write( "testmodel_global" )

print("use the interpolation of the global function")
u = grid2d.interpolate(velocityGlobal, "velocity", polorder=1)
m.setvelocity(u)
s.solve()
out_gl = grid2d.vtkOutput()
out_gl.add(s.solution())
out_gl.add(u)
out_gl.write( "testmodel_interpolation" )

print("use the discrete solution to some other scheme (vector valued)")
c0 = sympy.cos(2*math.pi*vecmodel.x0)
c1 = sympy.cos(2*math.pi*vecmodel.x1)
exact = [c0*c1, c0*c1]
u = vecmodel.trialFunction()
v = vecmodel.testFunction()
a = (ufl.inner(u,v) + ufl.inner(ufl.grad(u),ufl.grad(v)))*ufl.dx(0)
start_time = timeit.default_timer()
vecmodel.generateFromExact(a,exact)
VecModel = vecmodel.makeAndImport(grid2d)
print("Building VecModel took ", timeit.default_timer() - start_time, "s")
vecm = VecModel.get()
vecs = scheme.scheme("FemScheme", grid2d, vecm, "vector-valued", polorder=2)
vecs.solve()
sol=vecs.solution()
m.setvelocity(sol)
s.solve()
out_sol = grid2d.vtkOutput()
out_sol.add(s.solution())
out_sol.add(sol)
out_sol.write( "testmodel_df" )

print("use a 'local function adapter'")
class LocalExpr:
    def __init__(self,df):
        self.df = df
        self.dimR = 2
    def evaluate(self,en,x,r):
        # y = en.geometry().position(x)
        y = en.geometry().position( [1./3.,1./3.] )
        rr = self.df.evaluate(en,x)
        r[0] =  rr[0] * (1-y[1])**2
        r[1] =  rr[1] * (1-y[1])**2
localVelo = LocalExpr(sol)
velocityLocal = grid2d.getLocal( "local_velocity", localVelo )
m.setvelocity(velocityLocal)
s.solve()
out_locsol = grid2d.vtkOutput()
out_locsol.add(s.solution())
out_locsol.add(velocityLocal)
out_locsol.write( "testmodel_local" )
velocityLocal = "remove"

print("using a different 'local function adapter'")
class LocalExprA:
    def __init__(self,df):
        self.df = df
        self.dimR = 2
    def evaluate(self,en,x,r):
        y = en.geometry().position(x)
        # FieldMatrix not exported to python
        dx = self.df.jacobian(0,en,x)
        dy = self.df.jacobian(1,en,x)
        r[0] = -dy[0]
        r[1] =  dx[0]
localVelo = LocalExprA(sol)
# problem: when using velocityLocal = gf.getLocal( localVelo ) we get a seg fault
# when calling m.setVelocity(velocityLocal) - must be fixed
velocityLocalA = grid2d.getLocal( "nabla_x_u", localVelo )
m.setvelocity(velocityLocalA)
s.solve()
out_locsolA = grid2d.vtkOutput()
out_locsolA.add(s.solution())
out_locsolA.add(velocityLocalA)
out_locsolA.write( "testmodel_rotation" )

print("FINISHED")
