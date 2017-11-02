#!/bin/bash

mkdir /tmp/dune
cd /tmp/dune

export CMAKE_FLAGS=" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=TRUE \
  -DDUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS=TRUE \
  -DDUNE_PYTHON_INSTALL_LOCATION=system \
"
wget -qO - https://gitlab.dune-project.org/core/dune-common/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/core/dune-geometry/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/core/dune-grid/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/core/dune-istl/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/core/dune-localfunctions/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/extensions/dune-alugrid/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/extensions/dune-spgrid/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/staging/dune-uggrid/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/staging/dune-python/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/dune-fem/dune-fem/repository/archive.tar.gz?ref=releases/2.6 | tar xz
wget -qO - https://gitlab.dune-project.org/dune-fem/dune-fempy/repository/archive.tar.gz?ref=releases/2.6 | tar xz

patch -d dune-common-* -p 1 < $(dirname $0)/documentation.patch

./dune-common*/bin/dunecontrol all
./dune-common*/bin/dunecontrol make install

cd /
rm -rf /tmp/dune
