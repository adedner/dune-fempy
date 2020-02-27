from dune import create
from dune.fem.function import integrate
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, outer, dot,triangle, div, grad, dx, as_vector, transpose, Identity
from dune.ufl import NamedConstant, DirichletBC
from ufl import cos, sin, exp, sqrt
from math import pi,log10
from ufl import *

class AltComputeVel:
    # def __init__(self,spc,forcing,bcs,Re,velCorrectionStep=False,withCombSpace=False,useDG=False,withStab=False):
    def __init__(self,model,method,velCorrectionStep=False):
        spc      = model.spc
        Re       = model.Re
        forcing  = model.f
        #velocity bnd conditions
        # bcs      = model.bcs
        if model.useCombSPACE == True:
            spcU = spc.components[0]
            spcP = spc.components[1]
            dimension = spcU.grid.dimension
        else:
            dimension = spc.grid.dimension

        trial = TrialFunction(spc)
        test  = TestFunction(spc)
        u     = as_vector([trial[0],trial[1]])
        v     = as_vector([test[0],test[1]])
        p     = as_vector([trial[2]])
        q     = as_vector([test[2]])
        x     = SpatialCoordinate(spc)
        mu    = NamedConstant(spc, "mu")
        nu    = NamedConstant(spc, "nu")
        time  = NamedConstant(spc, "time")     # current time
        n = FacetNormal(spc.cell())
        # he = MaxFacetEdgeLength(uflSpace.cell())('+') # this is wrong
        # hT = MaxCellEdgeLength(uflSpace.cell())
        # hT = CellVolume(spcU.cell())
        # hF = FacetArea(spc.cell())
        # # he = FacetArea(uflSpace.cell()) / Min( avg('+'), avg('-') )
        # heInv = hF / avg( hT )
        #conditional(x[0]<-0.999,1.,0. )
        outflowVeldirich = conditional(x[0]>2.1999,0.,1. )

        dotdirich = conditional(x[0]<1e-3,1.,0. )
        pressdirich = conditional(x[0]>2.199,1.,0. )
        # dotdirich = conditional(x[0]+1<1e-3,1.,0. ) + conditional(x[1]+1<1e-3,1.,0. )
        # dotdirich =  conditional(x[0]<-0.999,1.,0.)

        dirichlet = dotdirich
        # dirichlet = conditional(x[0]+1<1e-3,1.,0. ) + conditional(x[0]>0.009 ,1.,0. )  + conditional(x[1]+1<1e-3 ,1.,0. )  + conditional(x[1]>0.009 ,1.,0. )

        # TODO remove gamma_t
        # gamma_t = exp(-2.*pi*pi*(1./Re)*time)
        #
        # u_d     = as_vector([1.5/(0.41*0.41)*(0.41-x[1])*4*x[1]*conditional(x[0]<1e-3,1.,0.),0] )
        # exact_u = as_vector([1.5/(0.41*0.41)*(0.41-x[1])*4*x[1]*conditional(x[0]<1e-3,1.,0.),0,0])


        Um = 1.5
        H  = 0.41

        u_d     = as_vector( [1.5/(0.41*0.41)*(0.41-x[1])*4*x[1],0] )

        # TODO remove gamma_t
        gamma_t = exp(-2.*pi*pi*(1./model.Re)*time)

        # u_d     = as_vector( [(1-x[1]*x[1])* conditional(x[0]<1e-3,1.,0.),0] )
        exact_u = as_vector([1.5/(0.41*0.41)*(0.41-x[1])*4*x[1]* conditional(x[0]<1e-3,1.,0.),0,0])

        self.oldSolution = spc.interpolate( [0,0,0], name="oldSol")
        # self.oldVelo = spcU.interpolate([0,]*spcU.dimRange,name="oldU")
        # self.oldSolution = spc.interpolate([0,]*3,name="oldSol")
        self.oldVelo  = as_vector( [self.oldSolution[0],self.oldSolution[1]])
        self.oldPress =  self.oldSolution[2]

        Du       = 0.5 *(grad(u)+grad(u).T)
        oldDu    = 0.5 *(grad(self.oldVelo)+grad(self.oldVelo).T)
        norm2_Du = (inner(Du, Du))
        # b = (pow(d**(1) + (abs_du)**0.5, (pnb-2)/1)*inner(du, grad(v)) ) * dx
        # gradModel    = -1*inner( p[0]*Identity(dimension), grad(v) ) * dx

        if forcing is not None:
            f = inner(forcing,v)*dx
        else:
            f = 0

        if velCorrectionStep:
            f += (dot(self.oldVelo,v) - (1/mu)*inner(v, grad(self.oldPress))) * dx
        else:
             f += (mu*dot(self.oldVelo,v) -
             inner(grad(self.oldVelo), outer(v, self.oldVelo))) * dx



        if velCorrectionStep:
            tentVelmodel = dot(u,v)  * dx
        else:
            tentVelmodel  = (mu*dot(u,v) + 1/model.Re*2*inner(Du, grad(v))) * dx
            # tentVelmodel += -1*inner(div(v),p[0]) * dx
            # #tentVelmodel += -1*inner( p[0]*Identity(dimension), grad(v) ) * dx
            # tentVelmodel -= 1*inner(div(u),q[0]) * dx
            # if withStab:
            #     tentVelmodel += 2*inner(div(u),q[0]) * dx
            #     # tentVelmodel += 1e-2* (inner(grad(p[0]), grad(q[0]))) * dx
            #     # tentVelmodel +=  1e2*(inner((p[0] - self.oldPress), q[0])) * dx
            #     tentVelmodel +=  1e0*(inner((p[0] - self.oldPress), q[0])) * dx


        if model.useCombSPACE == True:
            if model.useDG == False:
                bndpart       =  (100)  * inner( u, v) *outflowVeldirich* ds
                #bndpart       +=  (100) * pressdirich * inner( p[0], q[0]) * ds
                # bndpart      +=  (10)  * dirichlet * inner(u, v) * ds
                bndpart      -=  (100)  * dirichlet * inner( as_vector([exact_u[0],exact_u[1]]), v) * ds
                #bndpart       -=  (100) * pressdirich * inner( 1e-50*p[0], q[0]) * ds

                tentVelmodel  +=  bndpart

        # if model.useDG == True:
        #         pena = 1e4
        #         bndpart       =  nu/Re * pena * 1 * inner( u, v) *outflowVeldirich*  ds
        #         bndpart      +=  nu/Re * pena * pressdirich * inner( p[0], q[0]) * ds
        #         # bndpart      +=  (10)  * dirichlet * inner(u, v) * ds
        #         # bndpart      -=  (10)  * dirichlet * inner( as_vector([exact_u[0],exact_u[1]]), v) * ds
        #         tentVelmodel  +=  bndpart
        #         dgpart       =   -1*nu/Re*  (inner(outer(jump(u), n('+')), avg(grad(v))) + inner(avg(grad(u)), outer(jump(v), n('+')))) * dS
        #         dgpart      +=   nu/Re * pena * inner(jump(u), jump(v)) * dS
        #         dgpart      -=   1*nu/Re * (inner(outer(u_d, n), grad(v)) + inner(grad(u_d), outer(v, n))) * dirichlet * ds
        #         dgpart      -=   nu/Re * pena * inner(u_d, v) * dirichlet * ds
        #         dgpart      -=   1*nu/Re * (inner(outer(u, n), grad(v)) + inner(grad(u), outer(v, n))) * dirichlet * ds
        #         dgpart      -=   nu/Re * pena * inner(u, v) * dirichlet * ds
        #
        #         dgpart      +=   -1*(1-nu)/Re*  (inner(outer(jump(self.oldVelo), n('+')), avg(grad(v))) + inner(avg(grad(self.oldVelo)), outer(jump(v), n('+')))) * dS
        #         dgpart      +=   (1-nu)/Re * pena * inner(jump(self.oldVelo), jump(v)) * dS
        #         dgpart      -=   1*(1-nu)/Re * (inner(outer(u_d, n), grad(v)) + inner(grad(u_d), outer(v, n))) * dirichlet * ds
        #         dgpart      -=   (1-nu)/Re *pena * inner(u_d, v) * dirichlet * ds
        #         dgpart      +=   1*(1-nu)/Re * (inner(outer(u, n), grad(v)) + inner(grad(u), outer(v, n))) * dirichlet * ds
        #         dgpart      +=   (1-nu)/Re *pena * inner(u, v) * dirichlet * ds
        #         tentVelmodel +=   dgpart


        solverParameter={"fem.solver.newton.linabstol": 1e-11,
                         "fem.solver.newton.linreduction": 1e-11,
                         "fem.solver.newton.tolerance": 1e-10,
                         "fem.solver.newton.verbose": "true",
                         "fem.solver.newton.linear.verbose": "false"}


        if model.useCombSPACE == True:
            if model.useDG == True:
                self.scheme  = create.scheme("galerkin",[tentVelmodel==f],spc, solver="gmres", parameters = solverParameter)
            else:
                self.scheme = create.scheme("galerkin", tentVelmodel==f, spc, solver="gmres", parameters = solverParameter)
        else:
            if model.useDG == True:
                self.scheme  = create.scheme("galerkin",[tentVelmodel==f],spc)
            else:
                model = create.model("elliptic", spc.grid, tentVelmodel == f,
                        DirichletBC(spc, [None,None,None], 2),       # bottom/top
                        DirichletBC(spc, [None,None,None], 4),       # bottom/top
                        DirichletBC(spc, [1.5/(0.41*0.41)*(0.41-x[1])*4*x[1],0,None], 3),             # left
                        DirichletBC(spc, [0,0,None], 1))             # right
                # model = create.model("elliptic", spc.grid, tentVelmodel == f,
                #         DirichletBC(spc, [0,0,None], 2),       # bottom/top
                #         DirichletBC(spc, [None,None,0], 4),       # bottom/top
                #         DirichletBC(spc, [1-x[1]*x[1],0,None], 3),             # left
                #         DirichletBC(spc, [0,0,None], 1))             # right
                # model = create.model("elliptic", spc.grid, tentVelmodel == f, DirichletBC(spc,[exact_u[0],exact_u[1],None],1))
                self.scheme = create.scheme("h1", model, spc, parameters=solverParameter)


    def prepare(self,time=0,mu=0,nu=1,oldSolution=None,oldSolution1=None):
        # self.scheme.model.time = time
        self.scheme.model.mu = mu
        # self.scheme.model.nu = nu

        self.oldSolution.assign(oldSolution)
        self.oldVelo = as_vector( [self.oldSolution[0],self.oldSolution[1]])
        self.oldPress =  self.oldSolution[2]

    def solve(self,target):
        self.scheme.solve(target=target)
