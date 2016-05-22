#ifndef DUNE_FEMPY_PY_GRID_HH
#define DUNE_FEMPY_PY_GRID_HH

#include <string>
#include <utility>

#include <dune/fempy/py/grid/gridview.hh>
#include <dune/fempy/py/grid/hierarchical.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerGrid
    // ------------

    template< class GridView >
    void registerGrid ( pybind11::module module )
    {
      registerHierarchicalGrid< typename GridView::Grid >( module );
      module.def( "readDGF", &readDGF< typename GridView::Grid > );
      module.def( "makeSimplexGrid", &makeSimplexGrid< typename GridView::Grid > );

      registerGridView< GridView >( module, "LeafGrid" );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HH
