#ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
#define DUNE_FEMPY_PY_GRID_GRIDPART_HH

#include <string>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fempy/py/grid/range.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/py/grid/vtk.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerGridView
    // ----------------

    template< class GridView >
    pybind11::class_< GridView > registerGridView ( pybind11::handle scope, const char *name )
    {
      typedef typename GridView::Grid Grid;

      const int dim = GridView::dimension;

      pybind11::class_< GridView > cls( scope, name );
      cls.def( "__init__", [] ( GridView &instance, Grid &grid ) {
          new (&instance) GridView( grid );
        }, pybind11::keep_alive< 1, 2 >() );

      registerPyGridViewRange< GridView, 0 >( cls, "Elements" );
      cls.def_property_readonly( "elements", [] ( pybind11::object gridView ) {
          return PyGridViewRange< GridView, 0 >( gridView.template cast< const GridView & >(), gridView );
        } );

      registerPyGridViewRange< GridView, dim >( cls, "Vertices" );
      cls.def_property_readonly( "vertices", [] ( pybind11::object gridView ) {
          return PyGridViewRange< GridView, dim >( gridView.template cast< const GridView & >(), gridView );
        } );

      cls.def( "__repr__", [ name ] ( const GridView &gridView ) -> std::string {
          return std::string( name ) + " with " + std::to_string( gridView.indexSet().size( 0 ) ) + " elements";
        } );

      cls.def_property_readonly( "hierarchicalGrid", [] ( GridView &gridView ) -> const Grid & { return gridView.grid(); } );

      cls.def( "size", [] ( const GridView &gridView, int codim ) { return gridView.indexSet().size( codim ); } );

      registerVTKWriter< GridView >( cls );
      cls.def( "vtkWriter", [] ( const GridView &gridView ) {
          return new VTKWriter< GridView >( gridView );
        }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "globalGridFunction", defGlobalGridFunction< GridView >( cls, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ) );
      cls.def( "localGridFunction", defLocalGridFunction< GridView >( cls, "LocalGridFunction", std::make_integer_sequence< int, 11 >() ) );

      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
