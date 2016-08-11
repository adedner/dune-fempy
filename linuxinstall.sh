#!/bin/sh
if ! which git >/dev/null; then
    echo 'Git must be installed to run this script. Install with "sudo apt install git"'
    exit 1
fi
if ! which cmake >/dev/null; then
    echo 'Cmake must be installed to run this script. Install with "sudo apt install cmake"'
    exit 1
fi

echo 'This script is for installing dune-fempy on Linux

It will install the following python modules: 
pip setuptools mpi4py sympy ufl sphinx

It will install Open MPI

It will also install DUNE with the following components:
dune-common dune-geometry dune-grid dune-istl dune-localfunctions dune-alugrid dune-fem dune-fempy'

read -p "Continue (y/n)? " choice
case "$choice" in 
  y|Y|yes|Yes )
        read -p "Do you want the default python version to be set to python3.5m in .bashrc (recommended) (y/n)? " choice2
        case "$choice2" in
          y|Y|yes|Yes )
            alias python=/usr/bin/python3.5m
            echo 'alias python=/usr/bin/python3.5m' >> $HOME/.bashrc ;;
          n|N|no|No ) ;;
          * ) echo "Please choose y or n" 
              exit 1 ;;
        esac
        read -p "Please enter a folder name (default is dune): " folder_name
        if [ -z "$folder_name" ]; then
          echo 'Using default name "dune"'
          folder_name='dune'
        fi
        echo "Making directory: $folder_name"
        mkdir $HOME/$folder_name
        echo 'Installing python modules'
        wget https://bootstrap.pypa.io/get-pip.py && python3.5m get-pip.py --user
        rm get-pip.py
        export PATH=$PATH:$HOME/.local/bin
        echo 'export PATH=$PATH:'$HOME'/.local/bin' >> $HOME/.bashrc
        if ! which pip3.5 >/dev/null; then
            echo 'pip3.5 failed to install. Check whether python3.5 is available (or use alternative version)'
            exit 1
        fi
        pip3.5 install -U pip
        pip3.5 install --user setuptools
        pip3.5 install --user mpi4py
        pip3.5 install --user sphinx
        git clone https://bitbucket.org/fenics-project/ufl.git $HOME/.local/lib/python3.5/site-packages/ufl
        cd $HOME/.local/lib/python3.5/site-packages/ufl
        python3.5m setup.py install  --user
        cd $HOME/$folder_name/
        read -p "Install Open MPI (y/n)? " choice3
        case "$choice3" in
          y|Y|yes|Yes )
            echo 'Installing Open MPI'
            wget https://www.open-mpi.org/nightly/v1.6/openmpi-1.6.6rc1r31736.tar.bz2
            tar xvjf openmpi-1.6.6rc1r31736.tar.bz2
            rm openmpi-1.6.6rc1r31736.tar.bz2
            mv openmpi-1.6.6rc1r31736 openmpi
            cd openmpi
            ./configure --prefix=$HOME/$folder_name/openmpi --disable-dlopen
            make all install
            export PATH=$PATH:$HOME/$folder_name/openmpi/bin
            export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/$folder_name/openmpi/lib"
            echo 'export PATH=$PATH:'$HOME/$folder_name'/openmpi/bin' >> $HOME/.bashrc
            echo 'export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:'$HOME/$folder_name'/openmpi/lib"' >> $HOME/.bashrc;;
          n|N|no|No ) ;;
          * ) echo "Please choose y or n" 
              exit 1 ;;
        esac
        read -p "Install DUNE (make sure cmake is installed first) (y/n)? " choice4
        case "$choice4" in
          y|Y|yes|Yes )
            cd $HOME/$folder_name
            git clone -b releases/2.4 https://gitlab.dune-project.org/core/dune-common.git dune-common
            git clone -b releases/2.4 https://gitlab.dune-project.org/core/dune-geometry.git dune-geometry
            git clone -b releases/2.4 https://gitlab.dune-project.org/core/dune-grid.git dune-grid
            git clone -b releases/2.4 https://gitlab.dune-project.org/core/dune-istl.git dune-istl
            git clone -b releases/2.4 https://gitlab.dune-project.org/core/dune-localfunctions.git dune-localfunctions
            git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git dune-alugrid
            git clone -b releases/2.4 https://gitlab.dune-project.org/dune-fem/dune-fem.git dune-fem
            git clone ssh://git@gitlab.dune-project.org:22022/michael.sghaier/dune-corepy.git dune-corepy
            git clone ssh://git@gitlab.dune-project.org:22022/dune-fem/dune-fempy.git dune-fempy
            echo 'CMAKE_FLAGS=" -DCMAKE_CXX_FLAGS='\''$OPTFLAGS'\'' \
                  -DDUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS=TRUE \
                  -DALLOW_CXXFLAGS_OVERWRITE=ON \
                  -DENABLE_HEADERCHECK=ON \
                  -DCMAKE_DISABLE_FIND_PACKAGE_LATEX=TRUE \
                  -DCMAKE_DISABLE_DOCUMENTATION=TRUE \
                  -DBUILD_SHARED_LIBS=TRUE \
                  -DCMAKE_PREFIX_PATH='\''$MODPATH/PETSC-dev;$MODPATH/Zoltan;$MODPATH/petsc;'\'' \
                  -DBoost_NO_BOOST_CMAKE=TRUE \
                  -DBoost_NO_SYSTEM_PATHS=TRUE \
                  -DBOOST_ROOT:PATHNAME='\'''$HOME/$folder_name'/boost'\'' \
                  -DBoost_LIBRARY_DIRS:FILEPATH='\'''$HOME/$folder_name'/boost/stage/lib'\'' \
                  -DPYTHON_EXECUTABLE:FILEPATH=/usr/bin/python3.5m \
                  -DPYTHON_INCLUDE_DIR:PATH=/usr/include/python3.5m \
                  -DPYTHON_LIBRARY:FILEPATH=/usr/lib64/libpython3.5m.so \
                  " ' > config.opts
            echo 'Installing DUNE' 
            ./dune-common/bin/dunecontrol --opts=config.opts all
            export PYTHONPATH="'$HOME/$folder_name'/dune-fempy/build-cmake/python:$PYTHONPATH"
            echo 'export PYTHONPATH="'$HOME/$folder_name'/dune-fempy/build-cmake/python:$PYTHONPATH"' >> $HOME/.bashrc
            ;;
          n|N|no|No ) ;;
          * ) echo "Please choose y or n" 
              exit 1 ;;
        esac
        ;;
  n|N|no|No ) exit 1 ;;
  * ) echo "Please choose y or n" 
      exit 1 ;;
esac
