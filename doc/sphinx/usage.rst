.. _usage:

.. raw:: html

  <h1> Usage guide </h1>

################################
Introduction
################################

In the :ref:`tutorial <tutorial>`, we looked at the basic functions available to Dune-Fempy in the context of Laplace equation. Here we will explain in more detail how we set up the various parts of a numerical problem on the python side, and the tools we have at our disposal. 

Behind all of the interface methods we use, the philosophy is that they are set up in a very similar way to the Dune-Fem structure. This should hopefully make the underlying C++ code transparent and easier to understand for a user of the python code, and vice versa. For more information about the C++ code, see the :ref:`advanced topics <advanced>` section.

.. contents::

################################
Running the examples
################################

First of all, once Dune-Fempy has been installed, a good way to check whether everything works is to see whether the demos are working properly. The demos are located in the demo folder in the main directory, but in order to run them, you have to go into the build-cmake folder. For instance, to run the laplace.py demo, you would type

.. code-block:: bash

  cd build-cmake
  make
  cd demo
  python laplace.pyc

The make command is only necessary if any changes are made to the files.

################################
Setting up a computational grid
################################

In Dune-Fempy the **grid** (somewhat self-explanatorily) refers to the grid used in the numerical method. It contains information about the mesh file, the dimension, and the Dune type that the grid takes. Grids, much like other parts of the problem such as the space and the scheme, can be set up easily in python using the database found in python/database/grid. This allows the user to specify grids from various parts of Dune that they want to use (more details on the topic of databases can be found in :ref:`Database approach <database>`). An example of this in python is the following

.. code-block:: python

  grid = dune.fem.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)

###############################################
Setting up a space
###############################################

In Dune-Fempy the **space** refers to the function space used in our finite element method. The space can be set up in python in an identical way to the grid as follows

.. code-block:: python

  space = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

###############################################
Setting up a mathematical model using UFL
###############################################

In Dune-Fempy, the **model** refers to the part of the problem that contains the weak form of the PDE and its boundary conditions. UFL is used to express the PDE, and from this we can generate a Dune model file. The module generation is done in the file python/dune/models/elliptic.hh.

Let us consider an example of UFL used to represent the Laplace equation in 2D. i.e. we consider the following PDE in weak form

.. math::

  \int uv + \nabla u\cdot\nabla v  =  \int f v

Here we let the left hand side of the equation be the bilinear form :math:`a(u,v)` and the right hand side be the linear functional :math:`b(v)`. Using the 2D grid we defined earlier, we can express this in UFL as follows

.. code-block:: python

  grid = dune.fem.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)

  uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1)
  u = TrialFunction(uflSpace)
  v = TestFunction(uflSpace)
  x = SpatialCoordinate(uflSpace.cell())

  f = cos(2*math.pi*x[0])*cos(2*math.pi*x[1])

  a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
  b = f * v[0] * dx

Once the parts of the model have been declared using the above, a python model object can be generated using

.. code-block:: python

  model = dune.models.elliptic.importModel(grid, a == b).get()

Here ``a``, ``b`` are given above as the LHS and RHS parts of the PDE. The full implementation of this model is given in :ref:`usageexample`.

Boundary conditions
-------------------

Additionally, boundary conditions can also be added to the model using UFL. Any *natural* boundary conditions (e.g. Neumann or Robin) can be added to the weak form directly by using a surface integral ds (instead of dx). On the other hand, *essential* boundary conditions can be added optionally using the **dirichlet** argument as follows

.. code-block:: python

  g1 = [cos(x[0]), sin(x[0])]
  g2 = [x[1], 3]
  model = dune.models.elliptic.importModel(grid, a == b, dirichlet = {1:[g1], 2:[g2]}).get()

Here ``1:[g1]`` tells us that the function ``g1`` is set on the boundary assigned to ``1`` in the mesh file, and similarly ``2:[g2]`` sets boundary ``2`` to ``g2``. Multiple Dirichlet boundary conditions can be individually assigned to different boundaries in this way.

Coefficients
------------

Suppose we want to create a model with a function that can be set to different values without remaking the model each time. This has the advantage of saving time if we want to run the same model with slightly different parameters. Additionally this allows us to easily set a function to a solution previously computed in the code. We can do this using the **Coefficient** variable. Consider the following example (found in demo/afem.py)

.. code-block:: python

  uflSpace = UFLSpace(2, 1)
  u = TrialFunction(uflSpace)
  v = TestFunction(uflSpace)
  x = SpatialCoordinate(uflSpace.cell())
  bnd_u = Coefficient(uflSpace)

  def exact(x):
      phi = math.atan2(x[1], x[0])
      if x[1] < 0:
          phi += 2*math.pi
      return [(x.two_norm2**(90./cornerAngle)) * sin(180./cornerAngle*phi)]

  a = inner(grad(u), grad(v)) * dx

  model = importModel(grid, a == 0, dirichlet={1:[bnd_u]}, tempVars=False).get()
  model.setCoefficient(bnd_u.count(), grid.globalGridFunction("bnd", exact))

Here we declare ``bnd_u`` to be a Coefficient, and then set it to be assigned as a Dirichlet boundary condition as shown previously. Then after creating the model, we can set ``bnd_u`` using ``setCoefficient`` to be equal to the function ``exact``.

.. _dunemodel:

Stand-alone Dune model generation
---------------------------------

It is possible to just create a C++ model file using UFL code for use within the Dune-Fem-Howto framework without using any of the other python interface tools. The advantage of this is to forgo the complicated process of manually writing a model file with functions for the source, flux, linSource, linFlux and so on. This can be done quite easily in the following way.

1. Create a UFL model file in a similar way to above. For examples of exactly what is required, see the models folder for reference.
2. Run the generateModel script in the build-cmake/demos directory. For example, to generate a model file for the transport equation example, you would run

  .. code-block:: bash

    python generateModel.pyc ../../models/equation.py

  Optionally you can add -m or -t to the call to make a python module, or test it with a FEM scheme.
3. Use the generated model file in conjuction with your own Dune code to make a method. The file is outputted to build-cmake/python/dune/generated using the name given in the UFL file (e.g. TransportModel.hh in this case).

################################
Setting up a numerical scheme
################################

In Dune-Fempy, the **scheme** contains information about the method used to solve the PDE. Just as before, schemes can be set up in a similar way to grids and spaces using the database found in python/database/scheme. An example of this in python is the following

.. code-block:: python

  scheme = dune.fem.create.scheme("FemScheme", space, model, "scheme")

Here *space* and *model* must both be previously defined, as shown above.

.. _usageexample:

################################
A full example
################################

Here we give a complete example for a problem that uses all the above methods. Other such examples can be found in the demo directory.

>> more complicated example goes here 
