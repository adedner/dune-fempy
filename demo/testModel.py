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

from ufl import *

import dune.models.elliptic
import dune.ufl
import dune.common as common
import dune.fem as fem
import dune.fem.space as space
import dune.fem.scheme as scheme
import dune.fem.function as gf

grid = fem.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)
grid.hierarchicalGrid.globalRefine(3)

############################################################################
# build a transport model -eps Laplace(u) + velocity.grad(u) + gamma u = f
# with unknown velocity field
############################################################################
scalarSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1)
vectorSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), grid.dimWorld)
u    = TrialFunction(scalarSpace)
v    = TestFunction(scalarSpace)
x    = SpatialCoordinate(scalarSpace.cell())
velo = Coefficient(vectorSpace)

# diff = model.coefficient('diffusion')
diff = [1+0.1*x[0]*x[0]]

a = ( (dot(velo,grad(u[0]))+0.01*u[0])*v[0] + diff[0]*inner(grad(u[0]),grad(v[0])) ) * dx
L = 10./(1.+(x[0]*x[0]+x[1]*x[1])**4 )  *  v[0]*dx
########
start_time = timeit.default_timer()
# model.setCoefficient("diffusion",[1+0.1*model.x0*model.x0])
print("Building TransportModel took ", timeit.default_timer() - start_time, "s")
m = dune.models.elliptic.importModel(grid, a == L).get()

############################################################################
# construct a largrange scheme for the transport problem
############################################################################
sp = space.create( "Lagrange", grid, dimrange=m.dimRange, polorder=1 )
solution = sp.interpolate([0,]*m.dimRange) # , name="solution")
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
velocityGlobal = grid.function("global_velocity", 3, globalExpr=expr1)
m.setCoefficient(velo, velocityGlobal)
s.solve( None, solution )
vtk = grid.vtkWriter()
solution.addToVTKWriter(vtk, common.DataType.PointData)
velocityGlobal.addToVTKWriter(vtk, common.DataType.PointVector)
vtk.write("testmodel_global");
#################
print("use the interpolation of the global function")
u = grid.interpolate(velocityGlobal, space="Lagrange", name="interpolated_velocity", polorder=1)
m.setCoefficient(velo, u)
s.solve( target = solution )
vtk = grid.vtkWriter()
solution.addToVTKWriter(vtk, common.DataType.PointData)
u.addToVTKWriter(vtk, common.DataType.PointVector)
vtk.write("testmodel_interpolation");
#################
print("use the discrete solution to some other scheme (vector valued)")
u    = TrialFunction(vectorSpace)
v    = TestFunction(vectorSpace)
x    = SpatialCoordinate(vectorSpace.cell())
a = (inner(u,v) + inner(grad(u),grad(v)))*dx
c0 = cos(2*pi*x[0])
c1 = cos(2*pi*x[1])
L = ( c0*c1*v[0] + c0*c1*v[1] )*dx

start_time = timeit.default_timer()
vecm = dune.models.elliptic.importModel(grid, a == L).get()
vecsp = space.create( "Lagrange", grid, dimrange=vecm.dimRange, polorder=2 )
vecs = scheme.create( "FemScheme", vecsp, vecm, "computed_velocity" )
# the following is equivalent to:
#   veloh = vecsp.interpolate([0,]*vecm.dimRange, name="velocity")
#   vecs.solve(target = veloh)
veloh = vecs.solve()
m.setCoefficient(velo, veloh)
s.solve( None, solution )
vtk = grid.vtkWriter()
solution.addToVTKWriter(vtk, common.DataType.PointData)
veloh.addToVTKWriter(vtk, common.DataType.PointVector)
vtk.write("testmodel_df");
#################
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
velocityLocal = grid.function( "local_velocity", 3, localExpr=localVelo )
m.setCoefficient(velo, velocityLocal)
s.solve( None, solution )
vtk = grid.vtkWriter()
solution.addToVTKWriter(vtk, common.DataType.PointData)
velocityLocal.addToVTKWriter(vtk, common.DataType.PointVector)
vtk.write("testmodel_local");
velocityLocal = "remove"
####################
print("using a different 'local function adapter'")
class LocalExprA:
    def __init__(self,df):
        self.df = df
        self.dimR = 2
    def __call__(self,en,x):
        y = en.geometry.position(x)
        jac = self.df.localFunction(en).jacobian(x)
        return [ -jac[0][1],jac[0][0] ]
velocityLocalA = grid.function( "nabla_x_u", 3, localExpr=LocalExprA(veloh ))
s = scheme.create( "FemScheme", solution, m, "transport" ) # here solution is used if solve method does not specify target
m.setCoefficient(velo, velocityLocalA)
s.solve()
vtk = grid.vtkWriter()
solution.addToVTKWriter(vtk, common.DataType.PointData)
velocityLocalA.addToVTKWriter(vtk, common.DataType.PointVector)
vtk.write("testmodel_rotation");

print("FINISHED")
