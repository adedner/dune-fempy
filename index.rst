.. dune-fempy documentation master file, created by
   Andreas Dedner on Mon Mai 20 2019.

######################################
Welcome to dune-fempy's documentation!
######################################

************
Introduction
************

.. toctree::
   :maxdepth: 2

This module brings python scripting support to `Dune`_.
It serves three purposes:

1. High level program control for solving partial differential equations
   using classes from the `Dune`_ core and from `Dune-Fem`_.
   The unified form language `UFL`_ is used to describe the mathematical
   model, all realizations of the `Dune`_ grid interface can be used to
   work with the domain tesselation, and the finite element spaces,
   operator, and solvers provided by `Dune-Fem`_ for the descritizations
   and solving steps. All of this is available within to be used in python
   scripts or through jupyter notebooks.
2. Rapid prototyping of new methods or other parts of a simulation is easy
   since the interfaces provided are very similar to the `Dune`_ C++
   interface. This makes it easy to transfer a working prototype from
   python (easy to develop) to C++ (high efficiency). Small C++ code
   snippets can be easy called from python using just in time compilation.
3. Rapid prototyping of new implementations of `Dune`_ interfaces. So
   new implementations of the `Dune`_ grid interface can be easily
   tested. For `Dune-Fem`_ developers, new grid views, discrete function spaces, and
   scheme classes following the `Dune-Fem-Howto`_ concept can be added and tested.

.. _Dune: http://www.dune-project.org
.. _Dune-Fem-Howto: http://dune.mathematik.uni-freiburg.de/doc/html-howto/
.. _Dune-Fem: http://www.dune-project.org/fem/index.html
.. _UFL: http://fenicsproject.org/documentation/ufl/1.0-beta2/ufl.html


*************
This Document
*************

1. A simple scalar, non linear time dependent partial differential equation
   is used to describe the basic concepts. This leads through the steps
   required to set up the problem, solve the system of equations, and
   visualize the results. After this introduction we discuss how to use
   different solver backends (including for example `scipy`_ and `petsc`_).
   Finally we provide more detail on how to use the `Dune`_ grid interface,
   attach data to the grid entities and define general grid functions.
2. The examples in this section build in the general concepts described
   in the first part,
   introducing only a few new features. The :ref:`scripts`
   used for each of these and the following examples can be downloaded
   and are hopefully useful as starting point for new projects.
3. Local grid refinement and coarsening is a central feature of
   `Dune-Fem`_. Here we show how to use it for stationary and
   time dependent problems.
4. Finally other projects are presented some of them developed by
   the authors of this document, some contributed by other
   users. If you have used this package then we would like to
   hear about it and would ask you to contribute to this
   chapter. Doing so is quite easy (see :ref:`contributing` for details).

.. _scipy: https://www.scipy.org
.. _petsc: https://www.mcs.anl.gov/petsc/
.. _petsc4py: https://bitbucket.org/petsc/petsc4py

.. toctree::
   :maxdepth: 1
   :caption: Installation (out of date)

   installation


.. toctree::
   :maxdepth: 3
   :caption: A first example and general concepts

   gettingstarted


.. toctree::
   :caption: Some further examples
   :maxdepth: 3

   furtherexamples


.. toctree::
   :maxdepth: 3
   :caption: Grid Adaptivity

   adaptivity


.. toctree::
   :maxdepth: 3
   :caption: Moving Grids

   moving


.. toctree::
   :maxdepth: 3
   :caption: Additional Projects

   furtherprojects


.. toctree::
   :maxdepth: 3
   :caption: Information and Resources

   additional
.. API reference <api/modules>


******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
