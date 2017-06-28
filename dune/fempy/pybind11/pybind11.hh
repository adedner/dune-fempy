#ifndef DUNE_FEMPY_PYBIND11_PYBIND11_HH
#define DUNE_FEMPY_PYBIND11_PYBIND11_HH

#include <dune/fempy/pybind11/gridfunction.hh>

#include <dune/corepy/pybind11/complex.h>
#if HAVE_EIGEN
#include <dune/corepy/pybind11/eigen.h>
#endif
#include <dune/corepy/pybind11/extensions.h>
#include <dune/corepy/pybind11/numpy.h>
#include <dune/corepy/pybind11/pybind11.h>
#include <dune/corepy/pybind11/stl.h>

#endif // #ifndef DUNE_FEMPY_PYBIND11_PYBIND11_HH