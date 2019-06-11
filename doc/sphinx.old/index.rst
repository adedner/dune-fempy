.. dune-fempy documentation master file, created by
   sphinx-quickstart on Mon Apr  4 17:13:09 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dune-fempy's documentation!
======================================

############
Introduction
############

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

############
Installation
############

.. toctree::
   :maxdepth: 1

   installation

##################################
A few examples (Jupyter notebooks)
##################################

.. toctree::
   :maxdepth: 3

   laplace
   demos

#####################################
Some more details (Jupyter notebooks)
#####################################

.. toctree::
   :maxdepth: 3

   details
   duneufl

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Notebooks and other files
=========================


:download:`create <create.ipynb>`

:download:`laplace-dg <laplace-dg.ipynb>`

:download:`crystal <crystal.ipynb>`

:download:`laplace-dirichlet <laplace-dirichlet.ipynb>`

:download:`extending <extending.ipynb>`

:download:`laplace-intro <laplace-intro.ipynb>`

:download:`grid-construction <grid-construction.ipynb>`

:download:`laplace-la <laplace-la.ipynb>`

:download:`gridfunctions <gridfunctions.ipynb>`

:download:`mcf <mcf.ipynb>`

:download:`laplace-adaptive <laplace-adaptive.ipynb>`

:download:`parameters <parameters.ipynb>`

:download:`laplace-coefficients <laplace-coefficients.ipynb>`

:download:`spiral <spiral.ipynb>`

:download:`ufl-bindings <uflbindings.ipynb>`

:download:`parameter file example <parameter>`

:download:`unit cube grid file <unitcube-2d.dgf>`

:download:`sphere grid file <sphere.dgf>`

