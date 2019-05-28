.. _installation:

#######################
Using Docker or Vagrant
#######################

.. todo:: mention dune-fem-dev

###########
From Source
###########

************
Requirements
************

The following dependencies are needed for Dune-Fempy:

* At least C++11 compatible C++ compiler (e.g. gcc 5.3 or later)

* python (3.4 or later - possibly also works with 2.7 but not guaranteed)

  * mpi4py
  * ufl       (strongly recommended)
  * petsc4py  (optional)

* Dune (release 2.6 or later)

  * dune-common
  * dune-geometry
  * dune-grid
  * dune-istl
  * dune-localfunctions
  * dune-fem
  * dune-corepy
  * dune-alugrid  (strongly recommended)

******************************
Building the Dune Core modules
******************************

.. todo:: Mention available deb packages and perhaps link to other tutorials?

****************************
Setting up the python module
****************************

.. todo:: mention `dune-python/bin/setup-py.sh` and other options 

***************
Troubleshooting
***************

* The compiler version needs to be 5.3 or later. This can be checked in terminal with ::

  $ g++ --version

  If your version is out of date, you will need to upgrade your system to use Dune

* It is possible that the python version may be an issue. The script uses python3.5m. If during the Dune installation you get the error

  .. code-block:: none

    fatal error: pyconfig.h: No such file or directory

  This can probably be fixed by installing additional python3.5 libraries with e.g. ::

  $ sudo apt-get install libpython3.5-dev

  If python3.5 is not available on your system, you can simply change 3.5 for another appropriate version everywhere in the script (e.g. 3.4 or 2.7 (untested)). Otherwise, consider upgrading your system.

* One other problem is that a default version of Open MPI may already be installed. This will lead to errors where Dune appears to be looking in the wrong directory for Open MPI (e.g. usr/lib/openmpi instead of the home directory where the script installs it). This can be solved by running ::

  $ make uninstall

  in the original MPI install directory, followed by removing the folder. It will then be necessary to reinstall Open MPI and Dune. It may also be necessary to direct mpi4py to the new MPI installation. It is possible to check whether this is a problem by running python and trying out 

  .. code-block:: python

    from mpi4py import MPI

  If it comes up with an error, this can be fixed by installing mpi4py manually using the following commands ::

  $ git clone https://bitbucket.org/mpi4py/mpi4py.git
  $ cd mpi4py
  $ python setup.py build --mpicc=/path/to/openmpi/bin/mpicc
  $ python setup.py install --user

