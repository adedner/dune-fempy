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

* If during the installation of Dune, you get the error 

  .. code-block:: none

    cmake: command not found

  Then cmake needs to be installed with e.g. :: 

  $ sudo apt-get install cmake
  
* It is possible that the python version may be an issue. The script uses python3.5m since that is the latest version available at the time of writing. If during the Dune installation you get the error

  .. code-block:: none

    fatal error: pyconfig.h: No such file or directory
  
  This can probably be fixed by installing python3.5 with e.g. :: 

  $ sudo apt-get install libpython3.5-dev
  
* One other problem is that a default version of Open MPI may already be installed. This will lead to errors where Dune appears to be looking in the wrong directory for Open MPI (e.g. usr/lib/openmpi instead of the home directory where the script installs it). This can be solved by running ::

  $ make uninstall
  
  in the original MPI install directory, followed by removing the folder. It will then be necessary to reinstall Open MPI and Dune. Lastly, it may be necessary to direct mpi4py to the new MPI installation. It is possible to check whether this is a problem by running python and trying out 
  
  .. code-block:: python
  
    from mpi4py import MPI
  
  If it comes up with an error, this can be fixed by installing mpi4py manually using the following commands ::
  
  $ git clone https://bitbucket.org/mpi4py/mpi4py.git
  $ cd mpi4py
  $ python setup.py build --mpicc=/path/to/openmpi/bin/mpicc
  $ python setup.py install --user

#################################
Warwick CSC Desktops
#################################

For students and staff of the University of Warwick, one possibility is installation on a CSC desktop. For more details about getting an account, see `here <http://www2.warwick.ac.uk/fac/sci/csc/facilities/>`_.

To automatically install Dune-fempy on a CSC system, run the following in the main directory.

.. code-block:: bash

  $ sh cscinstall.sh
