.. _installation:

Installation
============

Here we provide a short guide for installing Dune-fempy on Linux and on CSC desktops.

#################################
Linux
#################################

The following dependencies are needed for Dune-fempy:

* C++11 or later (should be automatic)

* python2.7 or later (should be automatic)

  * setuptools 
  * mpi4py 
  * sympy 
  * ufl
  * sphinx

* Dune2.4 or later

  * dune-common 
  * dune-geometry 
  * dune-grid 
  * dune-istl 
  * dune-localfunctions 
  * dune-alugrid 
  * dune-fem 

* Open MPI

The above can all be automatically installed using the install script found in the main folder, i.e.

.. code-block:: bash

  $ sh linuxinstall.sh
  
**Troubleshooting**

If during the installation of Dune, you get the error 

.. code-block:: none

  cmake: command not found

Then cmake needs to be installed with e.g. :: 

  $ sudo apt-get install cmake
  
It is possible that the python version may be an issue. The script uses python3.5m since that is the latest version available at the time of writing. If during the Dune installation you get the error

.. code-block:: none

  fatal error: pyconfig.h: No such file or directory
  
This can probably be fixed with e.g. :: 

  $ sudo apt-get install libpython3.5-dev

#################################
Warwick CSC Desktops
#################################

For students and staff of the University of Warwick, one possibility is installation on a CSC desktop. For more details about getting an account, see `here <http://www2.warwick.ac.uk/fac/sci/csc/facilities/>`_.

To automatically install Dune-fempy on a CSC system, run the following in the main directory.

.. code-block:: bash

  $ sh cscinstall.sh
