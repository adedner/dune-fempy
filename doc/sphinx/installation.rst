.. _installation:

#######################
Using Docker or Vagrant
#######################

We provide files for building a Dune docker or vagrant environment on a host
machine which can be used to develop both C++ and Python code based on Dune.
In both cases the environment has to be build locally and can then be run from
any folder on the host machine, i.e., a folder containing some Dune module
or some Python script. Access rights are set so that the docker/vagrant user
has the same rights as the user who build the docker/vagrant image
in the directory from where it was started; consequently new files
generated during the session are modifiable on the host and vice versa.
The Linux distribution used is based on a Ubuntu image
and contains most programs needed for shell based code development.

Simply download this
:download:`script<https://gitlab.dune-project.org/dune-fem/dune-fem-dev/raw/master/rundocker.sh>`.
When executing this scripts the first time the docker image will be
downloaded which will take some time. After the download is completed
the working folder will be mounted into the docker container
under ``/host``. This way the script can be executed (without requiring an
additional image download) in any folder containing the Python project to
be worked on. In addition the scripts and notebooks discussed in documentation are
made available under ``/dunepy/DUNE/dune-fempy/docs``. The git repositories
of all required Dune modules are cloned under ``/dunepy/DUNE`` and a Python
virtual environment is setup. Additional Python packages can be easily
installed using ``pip install`` and additional Dune modules can be added
using ``git clone``. After adding a new Dune module in ``/dunepy/DUNE`` run
``updateAll`` to configure the new Dune module and update the Dune Python
package.

A detailed description of the *dune-fem docker development environment* is
given `here<https://gitlab.dune-project.org/dune-fem/dune-fem-dev>`.

*******
Vagrant
*******

Run `vagrant up` to get started and then `vagrant ssh` in your working
directory to go into the container.
To update the container run `vagrant provision` which will reexecute the
`bootstrap.sh` file.


###########
From Source
###########

.. note::
   We strongly encourage the use of a python virtual environment and the
   following instructions are written assuming that a virtual environment is
   activated.

************
Requirements
************

The following dependencies are needed for Dune-Fem python binding:

* At least C++11 compatible C++ compiler (e.g. gcc 5.3 or later)
* python (3.4 or later - possibly also works with 2.7 but not guaranteed)

  * mpi4py
  * numpy and scipy (strongly recommended)
  * matplotlib      (strongly recommended)
  * ufl             (strongly recommended)
  * petsc4py        (recommended)

* Required Dune modules (release 2.6 or later)

  * dune-common (https://gitlab.dune-project.org/core/dune-common.git)
  * dune-geometry (https://gitlab.dune-project.org/core/dune-geometry.git)
  * dune-grid (https://gitlab.dune-project.org/core/dune-grid.git)
  * dune-python (https://gitlab.dune-project.org/staging/dune-python.git)
  * dune-fem (https://gitlab.dune-project.org/dune-fem/dune-fem.git)

* Recommended Dune modules (releases 2.6 or later)

  * dune-istl (https://gitlab.dune-project.org/core/dune-istl.git)
  * dune-localfunctions (https://gitlab.dune-project.org/core/dune-localfunctions.git)
  * dune-alugrid  (https://gitlab.dune-project.org/extensions/dune-alugrid.git)

******************************
Building the Dune Core Modules
******************************

.. todo:: Mention available deb packages and perhaps link to other tutorials?

After cloning all the repositories simply run

.. code:: bash

   ./dune-common/bin/dunecontrol --opts=config.opts all

where :download:`config.opts<config.opts>` is an optional configuration
file containing for example flags for the `cmake` process using `CMAKE_FLAGS=`.

.. todo:: we need to mention `CMAKE_POSITION_INDEPENDENT_CODE=TRUE` or `BUILD_SHARED_LIBS`

********************************
Building the Dune Python Package
********************************

After the build process has terminated (hopefully successfully) run

.. code:: bash

   ./dune-python/bin/setup-dunepy.py --opts=config.opts install

and you should be ready to go. Test the installation by opening a Python
terminal and running

.. code:: python

   from dune.grid import structuredGrid
   grid = structuredGrid([0,0],[1,1],[10,10])
   grid.plot()

If you have everything set up correctly (and have `matplotlib`) you should
get a figure of a structured grid...

.. note::
   The first time you construct an object of a specific realization of one
   of the Dune interfaces (e.g. here a structured grid),
   the just in time compiler needs to be invoked. This can take quite some
   time - especially for grid realizations. This needs to be done only once
   so rerunning the above code a second time (even using other parameters
   in the `structuredGrid` function) should execute almost instantaniously.

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

