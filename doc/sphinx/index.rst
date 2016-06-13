.. dune-fempy documentation master file, created by
   sphinx-quickstart on Mon Apr  4 17:13:09 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dune-fempy's documentation!
======================================

#################################
Introduction
#################################

.. toctree::
   :maxdepth: 2

This module brings python scripting support to `Dune`_. 
It serves three purposes:

1. High level program control for solving partial differential equations
   using classes from the `Dune`_ core and from `Dune-Fem`_. This is
   described in detail in section :ref:`usageintro`.
2. Rapid prototyping of the model classes used in the `Dune-Fem-Howto`_.
   These classes provide the mathematical description of the partial
   differential equation to solve. Dune-Fempy uses `UFL`_ as input language
   and generates the corresponding model class. At this stage Dune-Fempy
   only supports a subset of the full `UFL`_ language. 
   A simple python script provides
   the stand-alone option to generate model classes from `UFL`_ input and
   to do simple unit testing of the generated classes. For users and
   developers of code based on the `Dune-Fem-Howto`_ this script provides a
   fast and easy approach for generating code for new mathematical models.
   Details can be found in section :ref:`modelintro`
   of this documentation. 
3. Rapid prototyping of new implementations of `Dune`_ interfaces. At the
   moment a new implementation of the `Dune`_ grid interface class can be
   tested. For `Dune-Fem`_ developers, new scheme classes following the
   `Dune-Fem-Howto`_ concept can be added and tested. More details on this
   aspect can be found in section :ref:`constructionintro`.

.. _Dune: http://www.dune-project.org
.. _Dune-Fem-Howto: http://dune.mathematik.uni-freiburg.de/doc/html-howto/
.. _Dune-Fem: http://www.dune-project.org/fem/index.html
.. _UFL: http://fenicsproject.org/documentation/ufl/1.0-beta2/ufl.html

#################################
Installation notes
#################################

.. toctree::
   :maxdepth: 2

   installation

#################################
Basic usage
#################################

.. toctree::
   :maxdepth: 1

   usage

#################################
PDE description
#################################

.. toctree::
   :maxdepth: 1

   model

#################################
On the fly module construction
#################################

.. toctree::
   :maxdepth: 1

   construction

###################################
`C++ Class documentation`_
###################################

.. _C++ Class documentation: file:../../doxygen/html/annotated.html

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

