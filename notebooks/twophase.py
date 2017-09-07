from __future__ import print_function
try:
    get_ipython().magic(u'matplotlib inline # can also use notebook or nbagg')
except:
    pass
import dune.common as common
from dune.grid import cartesianDomain, string2dgf
import dune.fem as fem
from dune.fem import parameter
from dune.fem.plotting import plotPointData as plot

import dune.create as create

parameter.append({"fem.verboserank": 0, "istl.preconditioning.method": "ilu-0", "istl.preconditioning.iterations": 1, "istl.preconditioning.relaxation": 1.2})

domain = cartesianDomain([0,0],[0.9,0.65],[30,20])
grid = create.view("adaptive", "ALUCube", domain, dimgrid=2)
hgrid = grid.hierarchicalGrid
order = 3
spc = create.space("dglagrange", grid, dimrange=2, order=order, storage="istl")
solution     = spc.interpolate([0,0], name="solution")
solution_old = spc.interpolate([0,0], name="solution_old")

import math
import ufl
from ufl import *

from dune.ufl import Space, DirichletBC

uflSpace = Space((grid.dimGrid, grid.dimWorld), 2)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
tau   = Constant(uflSpace.cell())

#### constants
g    = [0,]*grid.dimWorld ; g[grid.dimWorld-1] = -9.810 # [m/s^2]
g    = as_vector(g)
r_w  = 1000.  # [Kg/m^3]
mu_w = 1.e-3  # [Kg/m s]
r_n  = 1460.  # [Kg/m^3]
mu_n = 9.e-4  # [Kg/m s]

lensDomain = conditional(abs(x[1]-0.49)<0.03,1.,0.)*\
             conditional(abs(x[0]-0.45)<0.11,1.,0.)
def lens(a,b):
    return a*lensDomain + b*(1.-lensDomain)
Phi   = lens(0.39, 0.40)             # [-]
K     = lens(6.64*1e-14, 6.64*1e-11) # [m^2]
s_wr  = lens(0.10, 0.12)             # [-]
s_nr  = lens(0.00, 0.00)             # [-]
theta = lens(2.00, 2.70)             # [-]
pd    = lens(5000., 755.)            # [Pa]

#### initial conditions
p_w0 = (0.65-x[1])*9810.       # hydrostatic pressure
s_n0 = 0                       # fully saturated # boundary conditions
#### boundary (using exp to localize to part of boundary)
inflow = conditional(abs(x[0]-0.45)<0.06,1.,0.)*\
         conditional(abs(x[1]-0.65)<1e-8,1.,0.)
# this was used in the dune user meeting paper but I found -0.075 in some talk...
J_n  = -5.137*1e-5
J_w  = u[1]-u[1] # 0 fails due to code generation bug
dirichlet = conditional(abs(x[0])<1e-8,1.,0.) +\
            conditional(abs(x[0]-0.9)<1e-8,1.,0.)
p_wD = p_w0
s_nD = s_n0

s_n  = u[1]
s_w  = 1.-s_n

# Brooks Corey law (with regularizated cut off)
def Abs(a): return conditional(a>0,a,-a)
# def Abs(a): return ln(cosh(a*100.))/100. # a smoothed out |.|
def Max(a, b): return (a+b+Abs(a-b))/2.
def Min(a, b): return (a+b-Abs(a-b))/2.
# def cutOff(a): return Min(Max(a,0.005),0.995)
# compute a smooth cut off for over/under shoots
def cutOff(a): return Min(Max(a,0.0),1.0)

s_we = (s_w-s_wr)/(1.-s_wr-s_nr)
s_ne = (s_n-s_nr)/(1.-s_wr-s_nr)
s_we = cutOff(s_we)
s_ne = cutOff(s_ne)
kr_w = s_we**((2.+3.*theta)/theta)
kr_n = s_ne**2*(1.-(1.-s_ne)**((2.+theta)/theta))

p_c  = pd*s_we**(-1./theta)
p_w  = u[0]
p_n  = p_w + p_c

l_n  = kr_n / mu_n
l_w  = kr_w / mu_w
q_n  = 0
q_w  = 0

velocity_n = K*(grad(p_n)-r_n*g)
velocity_w = K*(grad(p_w)-r_w*g)

bulk_s1 = K*l_n*grad(p_c)
bulk_p1 = K*( (l_n+l_w)*grad(p_w) + l_n*grad(p_c) )
bulk_p2 = -K*( (r_n*l_n+r_w*l_w)*g )
bulk_p3 = q_w+q_n
bulk_s1 = K*l_n*grad(p_c)
bulk_s2 = K*(grad(p_w)-r_n*g)*l_n
bulk_s3 = q_n

#### bulk
form_p = ( inner(bulk_p1+bulk_p2,grad(v[0])) - bulk_p3*v[0] ) * dx
form_s = ( inner(bulk_s1+bulk_s2,grad(v[1])) - bulk_s3*v[1] ) * dx
#### boundary flux
form_p += (J_n+J_w) * v[0] * inflow * ds
form_s += J_n * v[1] * inflow * ds

######## dg terms
n, h = FacetNormal(uflSpace.cell()), MinFacetEdgeLength(uflSpace.cell())
# Q: shouldn't this involve at least K/mu_n so let's say we want 20 for
#    simple laplace then we should use a constant C satisfying
#    20 = C/(K/mu) = C 10^{11} 9.e-4 = C 10^8
#    so C = 2e-8 so much smaller then suggested by Birane, e.g.,
# penalty = 10. * order*(order+grid.dimWorld-1) / avg(h) * avg(K)/mu_n
# penalty_p_w = penalty
# penalty_s_n = penalty
penalty_p_w = 1.e-2 / avg(h)
penalty_s_n = 1.e-3 / avg(h)
#### skeleton
## penalty
form_p += penalty_p_w * inner(jump(p_w), jump(v[0])) * dS
form_s += penalty_s_n * inner(jump(s_n), jump(v[1])) * dS
## consistency
form_p -= inner(avg(bulk_p1+bulk_p2), outer(jump(v[0]), n('+'))) * dS
form_s -= inner(avg(bulk_s1+bulk_s2), outer(jump(v[1]), n('+'))) * dS
## symmetry
# form -= inner(outer(jump(p_w), n('+')), avg(K*(l_n+l_w)*grad(v[0]))) * dS
# form -= inner(outer(jump(s_n), n('+')), avg(K*l_n*grad(v[0]))) * dS
# form -= 1./avg(Phi) * inner(outer(jump(p_w), n('+')), avg(K*upwl_n*grad(v[1]))) * dS
# form -= 1./avg(Phi) * inner(outer(jump(s_n), n('+')), avg(K*upwl_n*grad(v[1]))) * dS

##### dirichlet conditions
## penalty
form_p += penalty_p_w * (p_w-p_wD) * v[0] * dirichlet * ds
form_s += penalty_s_n * (s_n-s_nD) * v[1] * dirichlet * ds
## consistency
form_p -= inner(bulk_p1+bulk_p2, outer(v[0], n)) * dirichlet * ds
form_s -= inner(bulk_s1+bulk_s2, outer(v[1], n)) * dirichlet * ds
## symmetry
# form -= inner(outer(p_w-p_wD, n), K*(l_n+l_w)*grad(v[0])) * dirichlet * ds
# form -= inner(outer(s_n-s_nD, n), K*l_n*grad(v[0])) * dirichlet * ds
# form -= inner(outer(p_w-p_wD, n), K*l_n*grad(v[1])) * dirichlet * ds
# form -= inner(outer(s_n-s_nD, n), K*l_n*grad(v[1])) * dirichlet * ds

form = Phi*(u[1]-solution_old[1])*v[1] * dx + form_p + \
        0.5*tau*(form_s + replace(form_s,{u:solution_old}))
model = create.model( "integrands", grid, form == 0)
model.setCoefficient(0,solution_old) # shouldn't be required - needs to be set during model generation

newtonParameter = {"linabstol": 1e-10, "linreduction": 1e-10, "tolerance": 5e-7,
        "verbose": "true", "linear.verbose": "false",
        "istl.gmres.restart": 50}
scheme = create.scheme("galerkin", spc, model,
         # solver=("suitesparse","umfpack"),
         parameters={"fem.solver.newton." + k: v for k, v in newtonParameter.items()})

maxLevel=3
tol=0.2
marker = common.Marker
def mark(element):
    solutionLocal = solution.localFunction(element)
    grad = solutionLocal.jacobian(element.geometry.domain.center)
    eta = grad[1].infinity_norm
    if eta > tol:
      return marker.refine if element.level < maxLevel else marker.keep
    else:
      return marker.coarsen

solution.interpolate(as_vector( [p_w0,s_n0] ))
dt = 1.
model.setConstant(tau,dt) # [s]  (5 in paper)
# the following is not really cool - one can't direcly output ufl expressions as yet
cutoff_h = create.function("ufl", grid, "cutoff", 5,
           as_vector([cutOff(solution[1]), cutOff(1.-solution[1])]))
velocity_nh = create.function("ufl", grid, "velocity_n", 5,
        replace(velocity_n,{u:solution}))
velocity_wh = create.function("ufl", grid, "velocity_w", 5,
        replace(velocity_w,{u:solution}))
print(cutoff_h.name)
count = 0
vtk = grid.writeVTK("twophaseB", pointvector=[velocity_nh,velocity_wh], pointdata=[solution, solution_old,cutoff_h], number=count)
if 0:
    for i in range(2):
        hgrid.mark(lambda e: marker.refine if e.geometry.center[1]>0.6 else marker.keep)
        fem.adapt(hgrid)
        fem.loadBalance(hgrid)
else:
    for i in range(2):
        print("pre adaptive (",i,"): ",grid.size(0),end="\n")
        solution_old.assign(solution)
        scheme.solve(target=solution)
        hgrid.mark(mark)
        fem.adapt(hgrid,[solution])
        fem.loadBalance(hgrid,[solution])
        solution.interpolate(as_vector( [p_w0,s_n0] ))

endTime = 2000.
t = 0.
++count
vtk.write("twophaseB", count) # why do we repeat the name here?

while t < endTime:
    print(t,grid.size(0),end="\n")
    solution_old.assign(solution)
    scheme.solve(target=solution)
    t += dt
    count += 1
    vtk.write("twophaseB", count) # why do we repeat the name here?
    hgrid.mark(mark)
    fem.adapt(hgrid,[solution])
    fem.loadBalance(hgrid,[solution])
