.. sectionauthor:: Andreas Dedner <a.s.dedner@warwick.ac.uk>, Robert Kl\ |oe|\ fkorn <robert.kloefkorn@iris.no>, Martin Nolte <nolte.mrtn@gmail.com>

##########
Grid Views
##########

When constructing a grid, the object returned to Python is always the so
called `LeafGridView`. Without any refinement this is simply a view on all
the elements of the grid. As soon as the grid is refined the leaf grid view changes
so that it always contains the `leaf` elements the grid, i.e., the elements
on the finest level. Since it is a read only view refinement is carried out
using the underlying `hierarchical grid`, i.e.,

.. code:: python

   grid.hierarchicalGrid.globalRefine(1)

For a given hierarchical grid one can use different views, i.e., a view on
all the elements of a given level:

.. include:: levelgridview_nb.rst

DUNE-FEM provides a number of additional views which will be discussed
further in this chapter:

*  `dune.fem.view.adaptiveLeafGridView`: this view should be used when
   the grid is supposed to be locally adapted. The view is still on the leaf
   elements of the grid but way data is attached to the entities of the grid
   is optimized for frequent grid changes as caused by local adaptivity. Its
   usage is shown in the following example.

*  `dune.fem.view.geometryGridView`: this is an example of a `meta` grid
   view which is constructed from an existing grid view and replaces some
   aspect - in this case the geometry of each element using a given grid
   function. This concept makes it easy to perform simulations for example on
   complex domains or on moving grids as shown in :ref:`Evolving Domains<geomGV>`.

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


The following two examples showcase adaptivity: the first one using a
residual a-posteriori estimator for an elliptic problem, the second one
shows adaptivity for a time dependent phase field model for crystal growth.
At the end of this section a dual weighted residual approach is used to
optimize the grid with respect to the error at a given point. While the
first two examples can be implemented completely using the available Python
bindings the final example requires using a small C++ snippet which is easy
to integrate into the Python code.

.. toctree::
   :maxdepth: 2

   laplace-adaptive_nb
   crystal_nb


################
Evolving Domains
################

As mentioned above DUNE-FEM provides a grid view that makes it easy to
exchange the geometry of each entity in the grid. To setup such a grid view
one first needs to construct a standard grid view, i.e., a `leafGridView`
and define a grid function over this view using for example a discrete
function, a UFL function, or one of the concepts described in the grid function section of
:ref:`General Concepts<concepts>` chapter. Note that the topology of the
grid does not change, i.e., how entities are connected with each other.
The following shows an example of how to change a grid of the unit square
into a grid of a diamond shape:

.. include:: geoview_nb.rst

By using a discrete function to construct a geometry grid view, it becomes
possible to simulate problems on evolving domains where the evolution is
itself the solution of the partial differential equation. We demonstrate
this approach based on the example of surface mean curvature flow first in
its simplest setting and then with the evolution of the surface depending
on values of a computed surface quantity satisfying a heat equation on the
surface:

.. todo:: add a coupled surface diffusion/evolution problem

.. toctree::
   :maxdepth: 2
   :name: geomGV

   mcf_nb

.. _algorithms:

#######################
Using C++ Code Snippets
#######################

.. todo:: link to `developers.rst`

In this section we demonstrate how it possible to use small piece of C++
code to either extend the existing functionality or improve the efficiency
of the code by moving code from Python to C++.  We have already seen how to
define grid functions using C++ code snippets in the
:ref:`Grid Function</concepts_nb.ipynb#Grid-Functions>` section of the
:ref:`General Concepts<concepts>` chapter.
In the following we will move parts of an algorithm from
Python to C++. The DUNE interfaces exported to
Python are very close to their C++ counterpart so that rapid prototyping of
new algorithms can be carried out using Python and then easily moved to
C++. This will be demonstrated in the following examples:

.. toctree::
   :maxdepth: 2

   mcf-algorithm_nb
   laplace-dwr_nb
   lineplot_nb
