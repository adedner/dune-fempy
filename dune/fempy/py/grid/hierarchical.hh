#ifndef DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH
#define DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH

#include <array>
#include <functional>
#include <list>
#include <map>
#include <memory>

#include <dune/grid/common/gridfactory.hh>

#include <dune/fempy/py/grid/entity.hh>
#include <dune/fempy/pybind11/functional.h>
#include <dune/fempy/pybind11/numpy.h>
#include <dune/fempy/pybind11/pybind11.h>
#include <dune/fempy/pybind11/stl.h>

namespace Dune
{

  namespace FemPy
  {

    // readDGF
    // -------

    template< class Grid >
    inline Grid *readDGF ( std::string dgf )
    {
      GridPtr< Grid > gridPtr( dgf );
      gridPtr->loadBalance();
      Grid *grid = gridPtr.release();
      return grid;
    }



    // makeSimplexGrid
    // ---------------

    template< class Grid, class float_t = typename Grid::ctype >
    inline Grid *makeSimplexGrid ( pybind11::array_t< float_t > points, pybind11::array_t< int > simplices )
    {
      typedef typename Grid::ctype ctype;

      GridFactory< Grid > factory;

      // insert points into factory

      pybind11::buffer_info bufPoints = points.request();
      if( (bufPoints.ndim != 2) || (bufPoints.shape[ 1 ] != Grid::dimensionworld) )
        throw std::invalid_argument( "points array must be of shape (*, " + std::to_string( Grid::dimensionworld ) + ")" );

      for( std::size_t i = 0; i < bufPoints.shape[ 0 ]; ++i )
      {
        const std::size_t offset = i * (bufPoints.strides[ 0 ] / sizeof( float_t ));
        FieldVector< ctype, Grid::dimensionworld > x;
        for( int j = 0; j < Grid::dimensionworld; ++j )
          x[ j ] = static_cast< ctype >( static_cast< float_t * >( bufPoints.ptr )[ offset + j * (bufPoints.strides[ 1 ] / sizeof( float_t )) ] );
        factory.insertVertex( x );
      }

      // insert simplices into factory

      pybind11::buffer_info bufSimplices = simplices.request();
      if( (bufSimplices.ndim != 2) || (bufSimplices.shape[ 1 ] != Grid::dimension+1) )
        throw std::invalid_argument( "simplices array must be of shape (*, " + std::to_string( Grid::dimension+1 ) + ")" );

      GeometryType type( GeometryType::simplex, Grid::dimension );
      std::vector< unsigned int > vertices( Grid::dimension+1 );
      for( std::size_t i = 0; i < bufSimplices.shape[ 0 ]; ++i )
      {
        const std::size_t offset = i * (bufSimplices.strides[ 0 ] / sizeof( int ));
        for( int j = 0; j <= Grid::dimension; ++j )
          vertices[ j ] = static_cast< int * >( bufSimplices.ptr )[ offset + j * (bufSimplices.strides[ 1 ] / sizeof( int )) ];
        factory.insertElement( type, vertices );
      }

      // create grid

      return factory.createGrid();
    }


    // registerHierarchicalGrid
    // ------------------------

    template< class Grid >
    inline auto registerHierarchicalGrid ( pybind11::handle scope )
    {
      registerGridEntities< Grid >( scope );

      pybind11::class_< Grid > cls( scope, "HierarchicalGrid" );
      cls.def( "__repr__", [] ( const Grid &grid ) -> std::string { return "HierarchicalGrid"; } );
      cls.def( "globalRefine", [] ( Grid &grid, int level ) { grid.globalRefine( level ); } );
      cls.def( "loadBalance", [] ( Grid &grid ) { grid.loadBalance(); } );

      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH
