# coding: utf-8

# # Adding new Interface Realizations and C++ Code [(Notebook)][1]
#
# [1]: _downloads/extending.ipynb

# In[ ]:

try:
    get_ipython().magic(u'matplotlib inline # can also use notebook or nbagg')
except:
    pass


# The aim of this section is to demonstrate how to extend the python package by either
# - providing a new realizations for an existing interface (e.g. a `space` or a `grid`)
# - implementing a new python function through some C++ code
#
# ### Adding realizations
# We start with a simple example of adding a discrete function space. The C++ implementation of this space is available in this module. We assume here that the

# In[ ]:

from dune.generator import hashIt
from dune.generator.generator import SimpleGenerator
from dune.fem.space import addAttr
def p1Bubble(gridview, dimrange, storage=None):
    # set the direct include path - all the include paths from the grid view used also need to be included
    includes = [ "dune/fem/space/p1bubble.hh" ] + gridview._includes
    dimw = gridview.dimWorld
    # new construct the actul C++ type name which is
    # template <class FunctionSpace, GridView >
    # Note: any exported C++ class contains (in addition to the include files needed)
    #       also the correct C++ type name of that class
    # The wrapper `Dune::FemPy::GridPart` is needed to make it possible to either use
    # `dune-fem` GridPart classes on Dune core `GridView`
    typeName = "Dune::Fem::BubbleElementSpace< " +      "Dune::Fem::FunctionSpace< double, double, " +          str(gridview.dimGrid) + ", " + str(dimrange) + " >, " +      "Dune::FemPy::GridPart< " + gridview._typeName + " > >"

    # Now add the information required for the binding - first the file containing the bindings
    includes = includes + ["dune/fempy/py/space.hh"]
    # now the module name (also used for the file)
    moduleName = "myspace_" + hashIt(typeName)
    generator = SimpleGenerator("Space", "Dune::FemPy")
    module = generator.load(includes, typeName, moduleName)
    addAttr(module, module.Space, "double", storage)
    return module.Space(gridview)


# In[ ]:

import dune.create as create
import dune.grid as grid
grid = create.grid("ALUConform", grid.cartesianDomain([0, 0], [1, 1], [4, 4]), dimgrid=2)
space = p1Bubble(grid,1)
print("size of the discrete space: ", space.size)
print("should be one per vertex and one per element: ", grid.size(0)+grid.size(2))


# In the `load` method it is also possible to provide additional constructors and methods which will be added to the generated python class. Assume for example that our space takes in addition to the `GridView` instance an additional integer:
# ~~~~
# constructor = ["[] ( " + typeName + " &self, "
#                        + gridview._typeName + " &gridview, int arg) {",
#                    "    new(&self) " + typename "
#                          ( gridPart< GridView >( gridView ), arg );"
#                    "  }, pybind11::keep_alive< 1, 2 >()"]
# ~~~~
# Then the new `load` call is simply
# ~~~~
# module = generator.load(includes, typeName, moduleName, constructor=constructor)
# ~~~~
# Note that the default constructor used for other spaces will not be added to the class if a custom constructor is passed in.
#
# Finally a short note on how to add the new space to the create mechnism. Let us assume that the method `p1Bubble` is part of a python package under the `dune` namespace package. For simplicity it is contained in the `__init__.py` file of the `dune.myspace` module.
# To create the space for example using
# ~~~
# dune.create.space("bubble",grid,1)
# ~~~
# add the following to the `__init__.py` file of the `dune.myspace` module:
# ~~~
# registry["space"] = {
#          "p1bubble"     : p1Bubble
#         }
# ~~~
# After installing `dune.myspace` e.g. using `pip` the new space can be constructed like any of the other spaces contained in the `dune.fem` package.

# ### Importing C++ code
# Now we discuss how to add stand alone C++ code using the `dune.generator` module.

# In[ ]:
