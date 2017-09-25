from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import inspect
import sys

from dune.common.compatibility import isString
from dune.generator.generator import SimpleGenerator
from dune.fem import function

from ._spaces import *

import ufl
import dune.ufl


def interpolate(space, func, name=None, **kwargs):
    """interpolate a function into a discrete function space

    Args:
        space: discrete function space to interpolate into
        func:  function to interpolate
        name:  name of the resulting discrete function

    Returns:
        DiscreteFunction: the constructed discrete function
    """
    if name is None:
        name = func.name
    return function.discreteFunction(space, name=name, expr=func, **kwargs)


def storageToSolver(storage):
    if storage == "adaptive" or storage == "fem":
        return "fem"
    elif storage == "istl":
        return "istl"
    elif storage == "numpy":
        return "numpy"
    elif storage == "eigen":
        return "eigen"
    elif storage == "petsc":
        return "petsc"

generator = SimpleGenerator("Space", "Dune::FemPy")

def addAttr(module, cls, field, storage):
    setattr(cls, "_module", module)
    setattr(cls, "field", field)
    setattr(cls, "interpolate", interpolate)
    setattr(cls, "numpyFunction", function.numpyFunction)

    if not storage:
        storage = str("fem")
    if isString(storage):
        import dune.create as create
        assert storageToSolver(storage), "wrong storage (" + storage + ") passed to space"
        storage = create.discretefunction(storageToSolver(storage))(cls)
    else:
        storage = storage(cls)
    setattr(cls, "storage", storage)

    from ufl.finiteelement import FiniteElementBase
    def uflSpace(self):
        space = dune.ufl.Space(self)
        return space
    cls.uflSpace = uflSpace
    def uflTrialFunction(self):
        trialFunction = ufl.TrialFunction(self.uflSpace())
        return trialFunction
    cls.uflTrialFunction = uflTrialFunction
    def uflTestFunction(self):
        testFunction = ufl.TestFunction(self.uflSpace())
        return testFunction
    cls.uflTestFunction  = uflTestFunction
    def uflSpatialCoordinate(self):
        spatialCoordinate = ufl.SpatialCoordinate(self.uflSpace().cell())
        # spatialCoordinate.duneSpace = self
        # find a way to get the space (or the grid) from the
        # spatialCoordinate?
        return spatialCoordinate
    cls.uflSpatialCoordinate = uflSpatialCoordinate
    def uflConstant(self, dimRange, name):
        if name:
            return dune.ufl.NamedConstant(self.uflSpace().cell(),dimRange,name)
        elif dimRange == 0:
            return ufl.Constant(self.uflSpace().cell())
        else:
            return ufl.VectorConstant(self.uflSpace().cell(), dim=dimRange)
    cls.uflConstant    = lambda self: uflConstant(self,0,None)
    cls.uflVectorConstant = lambda self,dimRange: uflConstant(self,dimRange,None)
    cls.uflNamedConstant  = lambda self, name, dimRange=0: uflConstant(self,dimRange,name)

fileBase = "femspace"

def module(field, storage, includes, typeName, *args):
    includes = includes + ["dune/fempy/py/space.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args)
    addAttr(module, module.Space, field, storage)
    return module
