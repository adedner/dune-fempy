"""
Functions for creating python modules and C++ classes for grids.

A small example - a more complete example can be found in testgrid.py

Examples:
    >>> from __future__ import print_function

    >>> # dunefempy modules
    >>> import pydunefem
    >>> print("import this module"); import grid
    import this module...

    >>> # just get the grid (only for testing - not used)
    >>> grid1 = dune.grid.leafGrid(str("../data/unitcube-1d.dgf"), str("OneDGrid") )

    >>> # get the full grid module and then the grid (module needed for grid # functions and output object)
    >>> yaspgrid2d = dune.grid.get(str("YaspGrid"), dimgrid=2)
    >>> grid2 = yaspgrid2d.LeafGrid(str("../data/unitcube-2d.dgf"))
"""

from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import sys

from types import ModuleType

from .generator import generator

class Generator(generator.Generator):
    def modifyIncludes(self, includes):
        return includes
    def modifyTypeName(self, typeName):
        return typeName+"::LeafGridView";

myGenerator = Generator("Grid")

def getGridType(grid, **parameters):
    """Return the grid type (using a function from database.py).
    """
    return myGenerator.getTypeName(grid, **parameters)

def addMethodsToLeafGrid(module,obj):
    setattr(obj, "_module", module)

def get(grid, **parameters):
    """Create a grid module using the grid-database.

    This function creates a python module called grid_xxx.py where xxx is a
    number derived from the grid type.
    It does this by fetching the typedef and includes from the grid-database
    using the given arguments.

    Example
        A generated grid module::

            #include <dune/grid/yaspgrid.hh>
            #include <dune/grid/io/file/dgfparser/dgfyasp.hh>

            typedef Dune::YaspGrid< 2, Dune::EquidistantCoordinates< double, 2 > > Grid;
            typedef Dune::Fem::AdaptiveLeafGridPart< Dune::YaspGrid< 2, Dune::EquidistantCoordinates< double, 2 > > > GridPart;

            BOOST_PYTHON_MODULE( gridc98ee9ab86a59240870d35a5c32ed4ab ) { registerDuneGrid(); }

    This would correspond to calling get("YaspGrid", dimgrid = 2). It also binds some functions to the created module.

    Args:
        grid (string): the identifier for the grid type to use
        parameters (kwargs): parameters used for fixing the grid type

            * dimgrid=int: dimension of the grid
            * dimworld=int: dimension of the world

    Returns:
        module: the newly created grid module

    """
    module = myGenerator.getModule(grid, **parameters)
    addMethodsToLeafGrid(module,module.LeafGrid)
    return module

def leafGrid(dgf, grid, **parameters):
    """Get a LeafGrid

    Do get() and create a C++ grid class (see grid.hh).

    Notes:
        This is equivalent to::

            gridmodule = get(grid,parameters)
            grid = gridmodule.LeafGrid(dgf)

    Args:
        dgf (string): the dgf file to use for construction the grid
        grid (string): the identifier for the grid type to use
        parameters (kwargs): parameters used for fixing the grid type

            * dimgrid=int: dimension of the grid
            * dimworld=int: dimension of the world

    Returns:
        LeafGrid: the constructed grid

    """
    if isinstance(grid, str):
        module = get(grid, **parameters)
    elif isinstance(grid, ModuleType):
        module = grid
    else:
        raise TypeError("leafGrid: 'grid' must be either a string or a module")

    ret = module.LeafGrid(module.readDGF(dgf))
    return ret

#############################################
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
