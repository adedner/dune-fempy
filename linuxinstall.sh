#!/bin/sh
if ! which git >/dev/null; then
    echo 'Git must be installed to run this script. Install with "sudo apt install git"'
    exit 1
fi
if ! which cmake >/dev/null; then
    echo 'Cmake must be installed to run this script.
    Install, e.g., with "sudo apt install cmake"'
    exit 1
fi

echo 'This script is for installing dune-fempy on Linux

It will install the following python modules:
pip setuptools mpi4py sympy ufl sphinx

It will install Open MPI

It will also install DUNE with the following components:
dune-common dune-geometry dune-grid dune-istl dune-localfunctions dune-alugrid dune-fem dune-corepy dune-fempy'

read -p "Continue (y/n)? " choice
case "$choice" in y|Y|yes|Yes )
        read -p "Please enter a folder name (default is dune): " folder_name
        if [ -z "$folder_name" ]; then
          echo 'Using default name "dune"'
          folder_name='dune'
        fi
        echo "Making directory: $folder_name"
        mkdir $HOME/$folder_name
        read -p "Install DUNE (make sure cmake is installed first) (y/n)? " choice4
        case "$choice4" in
          y|Y|yes|Yes )
            cd $HOME/$folder_name
            git clone -b releases/2.5 https://gitlab.dune-project.org/core/dune-common.git dune-common
            git clone -b releases/2.5 https://gitlab.dune-project.org/core/dune-geometry.git dune-geometry
            git clone -b releases/2.5 https://gitlab.dune-project.org/core/dune-grid.git dune-grid
            git clone -b releases/2.5 https://gitlab.dune-project.org/core/dune-istl.git dune-istl
            git clone -b releases/2.5 https://gitlab.dune-project.org/core/dune-localfunctions.git dune-localfunctions
            git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git dune-alugrid
            git clone -b releases/2.5 https://gitlab.dune-project.org/dune-fem/dune-fem.git dune-fem
            git clone https://gitlab.dune-project.org/staging/dune-corepy.git dune-corepy
            git clone ssh://git@gitlab.dune-project.org:22022/dune-fem/dune-fempy.git dune-fempy
            echo 'CMAKE_FLAGS=" -DCMAKE_CXX_FLAGS='\''$OPTFLAGS'\'' \
                  -DDUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS=TRUE \
                  -DALLOW_CXXFLAGS_OVERWRITE=ON \
                  -DENABLE_HEADERCHECK=ON \
                  -DCMAKE_DISABLE_FIND_PACKAGE_LATEX=TRUE \
                  -DCMAKE_DISABLE_DOCUMENTATION=TRUE \
                  -DBUILD_SHARED_LIBS=TRUE \
                  " ' > config.opts
            echo 'Installing DUNE'
            ./dune-common/bin/dunecontrol --opts=config.opts all
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
