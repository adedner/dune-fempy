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

Naturally we intend to solve this equation using a finite element method. Let us consider how this method is implemented in Dune-Fempy. First we want to set up the domain on a discretised mesh. We do this using the `leafGrid` method as follows 

.. code-block:: python

  grid = dune.fem.leafGrid(dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")

The first argument gives us the mesh used. In this case we set up a 16x16 square over :math:`[0,1] \times [0,1]`. The second argument, "ALUSimplexGrid", tells us the type of grid manager to use on the Dune side. The third argument sets the dimension of the grid which of course is 2, and finally `refinement="conforming"` is an optional argument that tells us whether our grid is conforming or non-conforming.

Next we want to set up the function space, which we do using the method `create.space`

.. code-block:: python

  spc = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

################################
Full example
################################

.. literalinclude:: ../../demo/laplace.py
   :language: python
