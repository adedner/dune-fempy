"""
  Documentation for the Grid Module
"""

from __future__ import print_function
# from functools import update_wrapper
import sys
import sympy,math
import ufl
import dune.models.femufl as duneuflmodel
import dune.fem.grid as grid
import dune.fem.scheme as scheme

yaspgrid = grid.get(str("YaspGrid"), dimgrid=2)
mygrid = grid.leafGrid(str("../../data/unitcube-2d.dgf"),"YaspGrid",dimgrid=2)

model    = duneuflmodel.DuneUFLModel(2,1,'Transport')
exact = [sympy.cos(2*math.pi*model.x0)*sympy.cos(2*math.pi*model.x1)]
u = model.trialFunction()
v = model.testFunction()
x = model.spatialCoordinate()
a = ( ufl.inner(ufl.grad(u[0]),ufl.grad(v[0])) ) * ufl.dx(0)
L = 10./(1.+(x[0]*x[0]+x[1]*x[1])**4 )  *  v[0]*ufl.dx(0)
model.generate(a,L,exact)
Model = model.makeAndImport(yaspgrid)
m = Model.get()
femscheme = scheme.scheme("FemScheme", mygrid, m, "transport", polorder=1)

##############################################################################

thismodule = sys.modules[__name__]

class Scheme( femscheme._module.Scheme ):
    """NEW
       The main leaf grid class"""

if __name__ == 'test':  # this is used by sphinx
    pass

if __name__ == '__main__': # protect against execution
    help(femscheme)        # das nicht - die neuen Methdoe auf der Instanz tauchen nicht auf
