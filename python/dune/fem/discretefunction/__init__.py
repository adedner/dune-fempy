from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib
import inspect
from dune.generator.generator import SimpleGenerator

from ._discretefunctions import *
from ._solvers import *

from dune.fem import function

try:
    import ufl
    from dune.ufl import GridFunction, expression2GF
except:
    pass

generator = SimpleGenerator("DiscreteFunction", "Dune::FemPy")

def addAttr(module, cls, storage):
    setattr(cls, "_module", module)
    setattr(cls, "_storage", storage)

fileBase = "femdiscretefunction"

def module(storage, includes, typeName, *args):
    includes = includes + ["dune/fempy/py/function/discrete.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args, bufferProtocol=True)
    addAttr(module, module.DiscreteFunction, storage)
    return module
