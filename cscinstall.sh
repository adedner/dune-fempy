#!/bin/sh
echo 'This script is for installing dune-fempy on Warwick CSC desktops

It will install the following python modules: 
pip setuptools mpi4py sympy ufl sphinx

It will also install DUNE with the following components:
dune-common dune-geometry dune-grid dune-istl dune-localfunctions dune-alugrid dune-fem dune-fempy'

read -p "Continue (y/n)? " choice
case "$choice" in 
  y|Y|yes|Yes )
        read -p "Do you want the default python version to be set to python3.4m in .bashrc (recommended) (y/n)? " choice2
        case "$choice2" in
          y|Y|yes|Yes )
	          alias python=/usr/bin/python3.4m
            echo 'alias python=/usr/bin/python3.4m' >> $HOME/.bashrc ;;
          n|N|no|No ) ;;
          * ) echo "Please choose y or n" 
              exit 1 ;;
        esac
        read -p "Configure Open MPI in .bashrc (y/n)? " choice3
        case "$choice3" in
          y|Y|yes|Yes )
            module load ompi/1.6.4/gnu/4.3.4
            export LD_LIBRARY_PATH="/warwick/openmpi/1.6.4/gnu/4.3.4/lib"
            echo 'module load ompi/1.6.4/gnu/4.3.4' >> $HOME/.bashrc
            echo 'export LD_LIBRARY_PATH="/warwick/openmpi/1.6.4/gnu/4.3.4/lib"' >> $HOME/.bashrc ;;
          n|N|no|No ) ;;
          * ) echo "Please choose y or n" 
              exit 1 ;;
        esac
        read -p "Please enter a folder name: " folder_name
        echo "Making directory: $folder_name"
        mkdir $HOME/$folder_name
        echo 'Installing python modules'
        wget https://bootstrap.pypa.io/get-pip.py && python3.4m get-pip.py --user
        rm get-pip.py
        export PATH=$PATH:$HOME/.local/bin
        echo 'export PATH=$PATH:'$HOME'/.local/bin' >> $HOME/.bashrc
        pip3 install -U pip
        pip3 install --user setuptools
        pip3 install --user mpi4py
        pip3 install --user sphinx
        git clone https://bitbucket.org/fenics-project/ufl.git $HOME/.local/lib/python3.4/site-packages/ufl
        cd $HOME/.local/lib/python3.4/site-packages/ufl
        python3.4m setup.py install  --user
        cd $HOME/$folder_name
        read -p "Install DUNE (y/n)? " choice4
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
            git clone ssh://git@gitlab.dune-project.org:22022/dune-fem/dune-fempy.git dune-fempy
            echo 'CMAKE_FLAGS=" -DCMAKE_CXX_COMPILER=g++-4.9 \
                  -DCMAKE_CXX_FLAGS='\''$OPTFLAGS'\'' \
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
                  -DPYTHON_EXECUTABLE:FILEPATH=/usr/bin/python3.4m \
                  -DPYTHON_INCLUDE_DIR:PATH=/usr/include/python3.4m \
                  -DPYTHON_LIBRARY:FILEPATH=/usr/lib64/libpython3.4m.so \
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
