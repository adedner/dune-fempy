############################################
Dynamic Local Grid Refinement and Coarsening
############################################

For refining and coarsening a grid locally the `dune.fem` module provides a
functon `adapt`. The storage of all discrete functions will be
automatically resized to accommodate the changes in the grid but the
resulting dof vector will not be initialized. To prolong and restrict data
from the old to the new grid, the corresponding discrete functions have to be
passed to the `dune.fem.adapt` method:

.. code:: python

    fem.adapt(u1,u2,...,uN)

The module `dune.fem` also provides a `globalRefine(level,*dfs)` method,
where a negative level globally coarsens the grid. If discrete functions
are passed in they will be prolong (restricted), the dof vectors of all
other dof vectors will be resized.

.. note::

   if the underlying storage of a discrete function is stored on the Python
   side as a numpy array, i.e., `vec = uh.as_numpy` was called, then access
   to `vec` will be undefined after a grid modification since the
   underlying buffer change will not have been registered.

The module `dune.fem` provides a function for marking elements for
refinement/coarsening:

.. code:: python

    def mark(indicator, refineTolerance, coarsenTolerance=0,
        minLevel=0, maxLevel=None):

where `indicator` is a grid function.
An element :math:`T` is marked for refinement if the value of ``indicator``
on :math:`T` is greater then ``refineTolerance`` and coarsened if the
value is less then ``coarsenTolerance``. The element :math:`T` is not
refined if its level is already at ``maxLevel`` and not coarsened if its
level it at ``minLevel``.  This method can for example be used to refine
the grid according to an equal distribution strategy by invoking

.. code:: python

    dune.fem.mark(indicator, theta/grid.size(0))

where `theta` is a given tolerance.

A layered Doerfler strategy is also available

.. code:: python

    def doerflerMark(indicator, theta, maxLevel=None, layered=0.05):


.. toctree::
   :maxdepth: 2

   laplace-adaptive_nb
   crystal_nb


############
Moving Grids
############

.. todo:: add some explanation on `GridParts`

.. toctree::
   :maxdepth: 2

   mcf_nb

.. _algorithms:

#######################
Using C++ Code Snipetts
#######################

.. todo:: add some explanation on `algorithms` and closeness of Python/C++ interface

.. literalinclude:: mcf-algorithm.py
   :pyobject: calcRadius

.. literalinclude:: radius.hh

.. toctree::
   :maxdepth: 2

   mcf-algorithm_nb
   lineplot_nb
