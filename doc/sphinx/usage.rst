.. _usageintro:

################################
Usage introduction
################################

As explained in the introduction, Dune-fempy provides a python interface for solving PDEs using Dune and Dune-Fem. Here we will explain how this can be used to set up the various parts of a numerical problem. For a complete example of how this works, see :ref:`usageexample`.

to do: explain how to run scripts (build-cmake -> make -> python test.pyc)

include a list of grids/spaces/schemes available in the database

perhaps a brief overview of the parts needed to set up and solve a scheme

.. contents::

################################
Setting up a computational grid
################################

In Dune-fempy the **grid** (somewhat self-explanatorily) refers to the grid used in the numerical method. It contains information about the mesh file, the dimension, and the Dune type that the grid takes. Grids, much like other parts of the problem such as the space and the scheme, can be set up easily in python using the database found in python/database/grid. This allows the user to specify grids from various parts of Dune that they want to use. An example of this in python is the following

.. code-block:: python

  grid = dune.fem.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)

An explanation of what the function leafGrid() does is given in the following docstring

.. autofunction:: dune.fem.grid.leafGrid()

Here the get() function does the following

.. autofunction:: dune.fem.grid.get()

###############################################
Setting up a space
###############################################

In Dune-fempy the **space** refers to the function space used in our finite element method. The space can be set up in python in an identical way to the grid as follows

.. code-block:: python

  space = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

This time the create.space() function is used, which does the following

.. autofunction:: dune.fem.create.space()

Here the get() function does the following

.. autofunction:: dune.fem.grid.get()

###############################################
Setting up a mathematical model using UFL
###############################################

In Dune-fempy, the **model** refers to the part of the problem that contains the weak form of the PDE and its boundary conditions. UFL is used to express the model, and from this we can generate a Dune model file. The module generation is done in the file python/dune/models/elliptic.hh. An example of UFL used to generate a model is given below.

.. code-block:: python

  uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1)
  u = TrialFunction(uflSpace)
  v = TestFunction(uflSpace)
  x = SpatialCoordinate(uflSpace.cell())

  f = cos(2*math.pi*x[0])*cos(2*math.pi*x[1])

  a = (inner(grad(u), grad(v)) + inner(u,v)) * dx(0)
  b = f * v[0] * dx(0)

Once the parts of the model have been declared using the above, a python model object can be generated using

.. code-block:: python

  model = dune.models.elliptic.importModel(grid, a == b).get()

Here *grid* must be previously defined in the code, and *a*, *b* are the LHS and RHS parts of the PDE given in UFL form. 

Additionally, boundary conditions can also be added to the model using UFL. Any natural boundary conditions (e.g. Neumann or Robin) can be added to the weak form directly by using a surface integral ds (instead of dx). On the other hand, *essential* boundary conditions can be added optionally as follows

.. code-block:: python

  g = [cos(x[0]), sin(x[0])]
  model = dune.models.elliptic.importModel(grid, a == b, dirichlet={1:g}).get()

Here the number 1 in {1:g} denotes which boundary it is assigned to in the mesh file. Dirichlet boundary conditions can be individually assigned to different boundaries in this way.

################################
Setting up a numerical scheme
################################

In Dune-fempy, the **scheme** contains information about the method used to solve the PDE. Just as before, schemes can be set up in a similar way to grids and spaces using the database found in python/database/scheme. An example of this in python is the following

.. code-block:: python

  scheme = dune.fem.create.scheme("FemScheme", space, model, "scheme")

Here *space* and *model* must both be previously defined, as shown above. An explanation of how scheme works is given in the following docstring

.. autofunction:: dune.fem.scheme.create()

Here as before, the function get() creates the C++ scheme class using the information given to it.

.. _usageexample:

################################
A full example
################################

Here we give a complete example for a problem that uses all the above methods. Other such examples can be found in the demo directory.

.. literalinclude:: ../../demo/laplace.py
   :language: python
