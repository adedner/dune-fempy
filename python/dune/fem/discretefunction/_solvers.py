from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

def solvers(includes, storage, operator):
    _, dfIncludes, dfTypeName, linearOperatorType, _ = storage
    includes += dfIncludes + ["dune/fempy/parameter.hh"]
    typeName = operator(dfTypeName,linearOperatorType)
    return includes, typeName


def femsolver(storage,solverType="gmres"):
    includes = ["dune/fem/solver/cginverseoperator.hh"]

    if solverType == "cg":
        operator = lambda df,_: "Dune::Fem::CGInverseOperator< " + df + " >"
    else:
        raise ValueError("wrong krylov solver - only cg available")

    includes, typeName = solvers(includes,storage,operator)
    return "fem",includes,typeName

def pardgsolver(storage,solverType="gmres"):
    includes = ["dune/fem/solver/pardginverseoperators.hh"]

    if solverType == "cg":
        operator = lambda df,_: "Dune::Fem::CGInverseOperator< " + df + " >"
    elif solverType == "gmres":
        operator = lambda df,_: "Dune::Fem::ParDGGeneralizedMinResInverseOperator< " + df + " >"
    elif solverType == "bicgstab":
        operator = lambda df,_: "Dune::Fem::ParDGBiCGStabInverseOperator< " + df + " >"
    else:
        raise ValueError("wrong krylov solver - only cg,gmres,bicgstab available")

    includes, typeName = solvers(includes,storage,operator)
    return "fem",includes,typeName

def oemfemsolver(storage,solverType="gmres"):
    includes = ["dune/fem/solver/oemsolver.hh"]

    if solverType == "cg":
        operator = lambda df,linop: "Dune::Fem::OEMCGOp<" + ",".join([df,linop]) + " >"
    elif solverType == "gmres":
        operator = lambda df,linop: "Dune::Fem::OEMGMRESOp< " + ",".join([df,linop]) + " >"
    elif solverType == "bicgstab":
        operator = lambda df,linop: "Dune::Fem::OEMBICGSTABOp< " + ",".join([df,linop]) + " >"
    elif solverType == "bicgsq":
        operator = lambda df,linop: "Dune::Fem::OEMBICGSQOp< " + ",".join([df,linop]) + " >"
    else:
        raise ValueError("wrong krylov solver - only cg, gmres, bicgstab, bicqsq available")

    includes, typeName = solvers(includes,storage,operator)
    return "oemfem",includes,typeName

def istlsolver(storage,solverType="gmres"):
    includes = ["dune/fem/solver/istlsolver.hh"]

    if solverType == "cg":
        operator = lambda df,linop: "Dune::Fem::ISTLCGOp<" + ",".join([df,linop]) + " >"
    elif solverType == "gmres":
        operator = lambda df,linop: "Dune::Fem::ISTLGMResOp< " + ",".join([df,linop]) + " >"
    elif solverType == "bicgstab":
        operator = lambda df,linop: "Dune::Fem::ISTLBICGSTABOp< " + ",".join([df,linop]) + " >"
    elif solverType == "minres":
        operator = lambda df,linop: "Dune::Fem::ISTLMINResOp< " + ",".join([df,linop]) + " >"
    else:
        raise ValueError("wrong krylov solver - only cg, gmres, bicgstab, minres available")

    includes, typeName = solvers(includes,storage,operator)
    return "istl",includes,typeName

def suitesparsesolver(storage,solverType="cg"):
    includes = ["dune/fem/solver/ldlsolver.hh", "dune/fem/solver/spqrsolver.hh", "dune/fem/solver/umfpacksolver.hh"]

    if solverType == "ldl":
        operator = lambda df,inop: "Dune::Fem::LDLOp<" + ",".join([df,linop]) + " >"
    elif solverType == "spqr_symmetric":
        operator = lambda df,linop: "Dune::Fem::SPQROp< " + ",".join([df,linop,"true"]) + " >"
    elif solverType == "spqr_nonsymmetric":
        operator = lambda df,linop: "Dune::Fem::SPQROp< " + ",".join([df,linop,"fakse"]) + " >"
    elif solverType == "umfpack":
        operator = lambda df,linop: "Dune::Fem::UMFPACKOp< " + ",".join([df,linop]) + " >"
    else:
        raise ValueError("wrong krylov solver - only ldl,spqr_symmetric,spqr_nonsymmetric,umfpack available")

    includes, typeName = solvers(includes,storage,operator)
    return "suitesparse",includes,typeName

def eigensolver(storage,solverType="bicgstab"):
    includes = ["dune/fem/solver/eigen.hh"]

    if solverType == "cg":
        operator = lambda df,_: "Dune::Fem::EigenCGInverseOperator< " + df + " >"
    elif solverType == "bicgstab":
        operator = lambda df,_: "Dune::Fem::EigenBiCGStabInverseOperator< " + df + " >"
    else:
        raise ValueError("wrong krylov solver - only cg,bicgstab available")

    includes, typeName = solvers(includes,storage,operator)
    return "eigen",includes,typeName
