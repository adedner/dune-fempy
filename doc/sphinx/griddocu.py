"""
  Documentation for the Grid Module
"""

from __future__ import print_function
# from functools import update_wrapper
import sys
import dune.fem.grid as grid

yaspgrid = grid.get(str("YaspGrid"), dimgrid=2)
mygrid = grid.leafGrid(str("../../data/unitcube-2d.dgf"),"YaspGrid",dimgrid=2)

thismodule = sys.modules[__name__]

def getLocal(name,func):
    pass
setattr(getLocal,"__doc__",grid.getLocal.__doc__)
def getGlobal(name,func):
    pass
setattr(getGlobal,"__doc__",grid.getGlobal.__doc__)

class LeafGrid( yaspgrid.LeafGrid ):
    """NEW
       The main leaf grid class"""
class VTKOutput( yaspgrid.VTKOutput ):
    """NEW
       Class for vtk output"""

if __name__ == 'test':  # this is used by sphinx
    pass

if __name__ == '__main__': # protect against execution
    help(mygrid)      # das nicht - die neuen Methdoe auf der Instanz tauchen nicht auf
