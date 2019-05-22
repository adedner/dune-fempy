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
   describe the domain tesselation, and the finite element spaces,
   operator, and solvers provided by `Dune-Fem`_ for the descritizations 
   and solving steps. All of this is available within to be used in python
   scripts or through jupyter notebooks.
2. Rapid prototyping of the model classes used in the `Dune-Fem-Howto`_.
   These classes provide the mathematical description of the partial
   differential equation to solve. Dune-Fempy uses `UFL`_ as input language
   and generates the corresponding model class.
   A simple python script provides
   the stand-alone option to generate model classes from `UFL`_ input and
   to do simple unit testing of the generated classes. For users and
   developers of code based on the `Dune-Fem-Howto`_ this script provides a
   fast and easy approach for generating code for new mathematical models.
   of this documentation.
3. Rapid prototyping of new implementations of `Dune`_ interfaces. So
   new implementations of the `Dune`_ grid interface can be easily
   tested. For `Dune-Fem`_ developers, new grid views, discrete function spaces, and
   scheme classes following the `Dune-Fem-Howto`_ concept can be added and tested.

.. _Dune: http://www.dune-project.org
.. _Dune-Fem-Howto: http://dune.mathematik.uni-freiburg.de/doc/html-howto/
.. _Dune-Fem: http://www.dune-project.org/fem/index.html
.. _UFL: http://fenicsproject.org/documentation/ufl/1.0-beta2/ufl.html

**************************
Installation (out of date)
**************************

.. toctree::
   :maxdepth: 1

   installation

************************************
A first example and general concepts
************************************

.. toctree::
   :maxdepth: 3

   gettingstarted

*********************
Some further examples
*********************

These examples build in the general concepts described above,
introducing only a few new features. The Python scripts
(jupyter notebooks) used for each example can be downloaded
and are hopefully useful as starting point for new projects.

.. toctree::
   :maxdepth: 3

   furtherexamples

***************
Grid Adaptivity
***************

Local grid refinement and coarsening is a central feature of
`Dune-Fem`_. Here we show how to use it for stationary and
time dependent problems.

.. toctree::
   :maxdepth: 3

   adaptivity

***************
Grid Adaptivity
***************

.. toctree::
   :maxdepth: 3

   moving

*******************
Additional Projects
*******************

Here are some other projects some of them developed by
the authors of this document, some contributed by other
users. If you have used this package then we would like to
hear about it and would ask you to contribute to this
chapter (see more details here).

.. toctree::
   :maxdepth: 3

   furtherprojects


Notebooks/Script and other files for download
=============================================

============================================= ============================================= =============================================
Example                                       Notebooks                                     Scripts
============================================= ============================================= =============================================
Bending beam (linear elasticity)              :download:`notebook <elasticity_nb.ipynb>`    :download:`script <elasticity.py>`
Crystal growth (phase field model)            :download:`notebook <crystal_nb.ipynb>`       :download:`script <crystal_nb.ipynb>`
============================================= ============================================= =============================================

:download:`unit cube grid file <unitcube-2d.dgf>`
:download:`sphere grid file <sphere.dgf>`


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


