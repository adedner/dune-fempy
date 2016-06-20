"""Compiles and solves different models (for code testing).
"""
############################################################################
# initialize
############################################################################
from __future__ import print_function
import sys,os,math, getopt
sys.path.append("../python")
import timeit
from mpi4py import MPI

import sympy
import ufl

import dune.models.femufl as duneuflmodel
import dune.fem.grid as grid
import dune.fem.space as space
import dune.fem.scheme as scheme
import dune.fem.function as gf

grid2d = grid.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)

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

dimR = m.dimRange
sp = space.create( "Lagrange", grid2d, dimrange=dimR, polorder=1 )
solution = sp.interpolate([0,]*dimR, name="velocity")
s = scheme.create( "FemScheme", sp, m, "transport" )

############################################################################
# Now use different way to solve this problem with given velocity field
############################################################################

print("define velocity field in global coordinates (using different approaches)")
expr1 = gf.MathExpression(["-(x1-0.5)","x0-1./2."])
#expr2 = gf.SympyExpression(["-(x1-0.5)","x0-1./2."])
def expr_func(x,r):
    r[0] = -(x[1]-0.5)
    r[1] = (x[0]-0.5)
#expr3 = gf.FuncExpression(2,expr_func)
velocityGlobal = grid2d.globalGridFunction("global_velocity", expr1)

m.setvelocity(velocityGlobal)
s.solve(solution)
vtk = grid2d.vtkWriter()
solution.addToVTKWriter(vtk, vtk.PointData)
velocityGlobal.addToVTKWriter(vtk, vtk.PointData)
vtk.write("testmodel_global");

print("use the interpolation of the global function")
u = grid2d.interpolate(velocityGlobal, space="Lagrange", name="velocity", polorder=1)
m.setvelocity(u)
s.solve(solution)
vtk = grid2d.vtkWriter()
solution.addToVTKWriter(vtk, vtk.PointData)
u.addToVTKWriter(vtk, vtk.PointData)
vtk.write("testmodel_interpolation");

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
vecsp = space.create( "Lagrange", grid2d, dimrange=vecm.dimRange, polorder=2 )
vecs = scheme.create( "FemScheme", vecsp, vecm, "vector-valued" )
veloh = vecsp.interpolate([0,]*vecm.dimRange, name="velocity")
vecs.solve(veloh)
m.setvelocity(veloh)
s.solve(solution)
vtk = grid2d.vtkWriter()
veloh.addToVTKWriter(vtk, vtk.PointData)
solution.addToVTKWriter(vtk, vtk.PointData)
vtk.write("testmodel_df");

print("use a 'local function adapter'")
class LocalExpr:
    def __init__(self,df):
        self.df = df
        self.dimR = 2
    def __call__(self,en,x):
        # y = en.geometry().position(x)
        y = en.geometry.position( [1./3.,1./3.] )
        rr = self.df.localFunction(en).evaluate(x)
        r = [rr[0] * (1-y[1])**2,\
             rr[1] * (1-y[1])**2]
        return r
localVelo = LocalExpr(veloh)
velocityLocal = grid2d.localGridFunction( "local_velocity", localVelo )
m.setvelocity(velocityLocal)
s.solve(solution)
vtk = grid2d.vtkWriter()
solution.addToVTKWriter(vtk, vtk.PointData)
velocityLocal.addToVTKWriter(vtk, vtk.PointData)
vtk.write("testmodel_local");
velocityLocal = "remove"

print("using a different 'local function adapter'")
class LocalExprA:
    def __init__(self,df):
        self.df = df
        self.dimR = 2
    def __call__(self,en,x):
        y = en.geometry.position(x)
        jac = self.df.localFunction(en).jacobian(x)
        return [ -jac[0][1],jac[0][0] ]
localVelo = LocalExprA(veloh)
# problem: when using velocityLocal = gf.getLocal( localVelo ) we get a seg fault
# when calling m.setVelocity(velocityLocal) - must be fixed
velocityLocalA = grid2d.localGridFunction( "nabla_x_u", localVelo )
m.setvelocity(velocityLocalA)
s.solve(solution)
vtk = grid2d.vtkWriter()
solution.addToVTKWriter(vtk, vtk.PointData)
velocityLocalA.addToVTKWriter(vtk, vtk.PointData)
vtk.write("testmodel_rotation");

print("FINISHED")
