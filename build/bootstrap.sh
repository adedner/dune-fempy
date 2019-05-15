#!/usr/bin/env bash

# This script installs all necessary packages to build DUNE-FemPy.

#Make sure that script exits on failure, and that all commands are printed
set -e
set -x

export DEBIAN_FRONTEND=noninteractive
# Make sure we have updated URLs to packages etc.
apt-get -y update

# system packages
# apt-get install -y ubuntu-desktop
apt-get install -y apt-utils
apt-get install -y ca-certificates
apt-get install -y pkg-config patch git vim gnuplot
apt-get install -y unzip gdb gzip tar wget time sudo
apt-get install -y build-essential gfortran gcc g++ libmpich-dev cmake cmake-curses-gui
apt-get install -y python3-venv python3-dev python3-setuptools python3-pip python3-wheel \
  python3-appdirs python3-configobj python3-requests python3-pandocfilters python3-tk
apt-get install -y libsuitesparse-dev petsc-dev
