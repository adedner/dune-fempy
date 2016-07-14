.. _tutorial:

.. raw:: html

  <h1> Tutorial </h1>

################################
Introduction
################################

As explained in the introduction, Dune-Fempy provides a python interface for solving PDEs using Dune and Dune-Fem. In this section we will go step-by-step through a simple example (the Laplace equation) to give an idea of what this interface looks like. For a more in-depth look at how the parts of Dune-Fempy work, see the :ref:`usage guide <usage>`.

.. contents::

################################
Laplace equation
################################

Let us consider the Laplace equation (to be more precise the Helmholtz equation) in a 2D domain :math:`\Omega` with boundary :math:`\Gamma`

.. math::

  \begin{gather}
  - \Delta u + u = f \quad \text{in } \Omega \\
  \nabla u \cdot \textbf{n} = 0 \quad \text{on } \Gamma
  \end{gather}

The weak form of this equation for a test function :math:`v`, is

.. math::

  \int_{\Omega} uv + \nabla u\cdot\nabla v \ dx =  \int_{\Omega} f v \ dx

Naturally we intend to solve this equation using a finite element method. Let us consider how this method is implemented in Dune-Fempy. First we want to set up the domain on a discretised mesh. We do this using the `leafGrid` method as follows.

.. code-block:: python

  grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", , refinement="conforming")

The first argument gives us the mesh used. In this case we set up a 16x16 square over :math:`[0,1] \times [0,1]` (though we can also give it a mesh file). The second argument, `"ALUSimplexGrid"`, tells Dune-Fempy the type of grid manager to use on the Dune side. The last two arguments are keyword arguments (note the = sign), hence the order they are given does not matter. `dimgrid=2` and `refinement="conforming"` tell us the dimension of the grid, and whether it is conforming or non-conforming.

Next we want to set up the function space, which we do using the method `create.space`.

.. code-block:: python

  spc = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

As before, the argument in quotes, `"Lagrange"`, tells us what Dune type to use when constructing the space. The second argument is a python grid module, which we defined above. `dimrange=1` and `polorder=2` give us the dimension of the range and the order of the finite elements we use.

The PDE itself can be expressed using UFL. In the above equation, we let the left hand side be the bilinear form :math:`a(u,v)` and the right hand side be the linear functional :math:`b(v)`. Then the code looks as follows.

.. code-block:: python

  uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1)
  u = TrialFunction(uflSpace)
  v = TestFunction(uflSpace)
  x = SpatialCoordinate(uflSpace.cell())

  f = cos(2*math.pi*x[0])*cos(2*math.pi*x[1])

  a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
  b = f * v[0] * dx

Once the parts of the model have been declared using the above, a python model object can be generated using `importModel`.

.. code-block:: python

  model = dune.models.elliptic.importModel(grid, a == b).get()

Note that we input the grid as an argument, and for the model itself we set the LHS and RHS to be equal.

Lastly it remains to define the method we use to solve the PDE. We do this by setting up the scheme as follows.

.. code-block:: python

  scheme = dune.fem.create.scheme("FemScheme", spc, model, "scheme")

Once again `FemScheme` is the Dunetype, `space` and `model` are previously defined, and "scheme" is the name we attach to the scheme.

Finally we can solve the model and output the data. To do the first we call solve on the scheme.

.. code-block:: python

  solution = scheme.solve()

And we create an output file using the method `writeVTK` on the grid.

.. code-block:: python

  grid.writeVTK("laplace", pointdata=[solution])

################################
Running the example
################################

The full python code for this example is given below.

.. literalinclude:: ../../demo/laplace.py
   :language: python

This example and others are located in the demo folder in the main directory. In order to run them, you have to go into the build-cmake folder and run the .pyc file. For instance, for this example you would type

.. code-block:: bash

  cd build-cmake
  make
  cd demo
  python laplace.pyc

Note that the make command is only necessary if any changes are made to the files.
