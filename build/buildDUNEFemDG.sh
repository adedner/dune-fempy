#!/bin/bash

# create necessary python virtual environment
if ! test -d $HOME/dune-env ; then
  python3.6 -m venv $HOME/dune-env
  source $HOME/dune-env/bin/activate
  pip install --upgrade pip
  pip install ufl numpy matplotlib scipy jupyter ipython mpi4py
echo "\
source $HOME/dune-env/bin/activate
export DUNE_LOG_FORMAT='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
export DUNE_LOG_LEVEL=CRITICAL
 " >> .bashrc
else
  source $HOME/dune-env/bin/activate
fi

#change appropriately, i.e. 2.6 or empty which refers to master
DUNEVERSION=

FLAGS="-O3 -DNDEBUG -funroll-loops -finline-functions -Wall -ftree-vectorize -fno-stack-protector -mtune=native"

DUNECOREMODULES="dune-common dune-istl dune-geometry dune-grid dune-localfunctions"
DUNEEXTMODULES="dune-python dune-alugrid"
DUNEFEMMODULES="dune-fem dune-fempy dune-fem-dg dune-vem"

#PY_CXXFLAGS=`python3-config --includes`
FLAGS="$FLAGS $PY_CXXFLAGS"

#PY_LDFLAGS=`python3-config --ldflags`

if ! test -d DUNE ; then
# we need to put the dune module into a subdirectory otherwise dune-py in
# the virtual env will be picked up during build of dune
mkdir DUNE
cd DUNE
# build flags for all DUNE and OPM modules
# change according to your needs
echo "\
DUNEPATH=`pwd`
BUILDDIR=build-cmake
USE_CMAKE=yes
MAKE_FLAGS=-j4
CMAKE_FLAGS=\"-DCMAKE_CXX_FLAGS=\\\"$FLAGS\\\"  \\
 -DDUNE_PYTHON_INSTALL_EDITABLE=TRUE \\
 -DADDITIONAL_PIP_PARAMS="-upgrade" \\
 -DCMAKE_LD_FLAGS=\\\"$PY_LDFLAGS\\\" \\
 -DCMAKE_POSITION_INDEPENDENT_CODE=TRUE \\
 -DDISABLE_DOCUMENTATION=TRUE \\
 -DCMAKE_DISABLE_FIND_PACKAGE_LATEX=TRUE\" " > config.opts

DUNEBRANCH=
if [ "$DUNEVERSION" != "" ] ; then
  DUNEBRANCH="-b releases/$DUNEVERSION"
fi

# get all dune modules necessary
for MOD in $DUNECOREMODULES ; do
  git clone $DUNEBRANCH https://gitlab.dune-project.org/core/$MOD.git
done

# get all dune extension modules necessary
for MOD in $DUNEEXTMODULES ; do
  if [ "$MOD" == "dune-alugrid" ]; then
    git clone $DUNEBRANCH https://gitlab.dune-project.org/extensions/$MOD.git
  else
    git clone $DUNEBRANCH https://gitlab.dune-project.org/staging/$MOD.git
  fi
done

# get all dune extension modules necessary
for MOD in $DUNEFEMMODULES ; do
  if [ "$MOD" == "dune-fem-dg" ]; then
    git clone -b advection-diffusion https://gitlab.dune-project.org/dune-fem/$MOD.git
  else
    git clone $DUNEBRANCH https://gitlab.dune-project.org/dune-fem/$MOD.git
  fi
done
else
cd DUNE
./dune-common/bin/dunecontrol git pull
rm -rf */build-cmake
fi

# build all DUNE modules using dune-control
./dune-common/bin/dunecontrol --opts=config.opts all

# install all python modules in the pip environment
./dune-python/bin/setup-dunepy.py --opts=config.opts install

##################################################################

cd ..
