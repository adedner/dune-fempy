#ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
#define DUNE_FEMPY_PY_GRID_FUNCTION_HH

#include <functional>
#include <string>
#include <tuple>
#include <utility>

#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fempy/py/grid/vtk.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerLocalFunction
    // ---------------------

    template< class LocalFunction >
    pybind11::class_< LocalFunction > registerLocalFunction ( pybind11::handle scope, const char *clsName = "LocalFunction" )
    {
      typedef typename LocalFunction::LocalCoordinateType LocalCoordinate;

      pybind11::class_< LocalFunction > cls( scope, clsName );

      cls.def_property_readonly( "dimRange", [] ( LocalFunction & ) -> int { return LocalFunction::RangeType::dimension; } );
      cls.def( "evaluate", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::RangeType value;
          lf.evaluate( x, value );
          return value;
        } );
      cls.def( "jacobian", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::JacobianRangeType jacobian;
          lf.jacobian( x, jacobian );
          return jacobian;
        } );

      return cls;
    }



    namespace detail
    {

      // registerGridFunction
      // --------------------

      template< class GridFunction >
      pybind11::class_< GridFunction > registerGridFunction ( pybind11::handle scope, const char *clsName = "GridFunction" )
      {
        typedef typename GridFunction::LocalFunctionType LocalFunction;
        typedef typename LocalFunction::EntityType Entity;

        pybind11::class_< GridFunction > cls( scope, clsName );

        registerLocalFunction< LocalFunction >( cls );

        cls.def( "__repr__", [] ( GridFunction &gf ) -> std::string {
            return "GridFunction< " + std::to_string( GridFunction::RangeType::dimension ) + " >(name = " + gf.name() + ")";
          } );

        cls.def_property_readonly( "dimRange", [] ( GridFunction &gf ) -> int { return GridFunction::RangeType::dimension; } );

        cls.def_property_readonly( "name", [] ( GridFunction &gf ) -> std::string { return gf.name(); } );
        cls.def_property_readonly( "grid", &GridFunction::gridView );

        cls.def( "localFunction", [] ( const GridFunction &gf, const Entity &entity ) -> LocalFunction {
            return gf.localFunction( entity );
          }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >() );

        cls.def( "addToVTKWriter", &addToVTKWriter< GridFunction >, pybind11::keep_alive< 1, 2 >() );

        return cls;
      }



    } // namespace detail



    // registerGridFunction
    // --------------------

    template< class GridFunction >
    pybind11::class_< GridFunction > registerGridFunction ( pybind11::handle scope, const char *clsName = "GridFunction" )
    {
      auto cls = detail::registerGridFunction< GridFunction >( scope, clsName );

      return cls;
    }




    namespace detail
    {

      // makePyGlobalGridFunction
      // ------------------------

      template< class GridView, int dimRange >
      auto makePyGlobalGridFunction ( const GridView &gridView, std::string name, pybind11::function evaluate, std::integral_constant< int, dimRange > )
      {
        typedef typename GridView::template Codim< 0 >::Geometry::GlobalCoordinate Coordinate;
        return simpleGridFunction( std::move( name ), gridView, [ evaluate ] ( const Coordinate &x ) {
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( x ) );
            return v.template cast< FieldVector< double, dimRange > >();
          } );
      }



      // registerPyGlobalGridFunction
      // ----------------------------

      template< class GridView, int dimRange >
      auto registerPyGlobalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyGlobalGridFunction( std::declval< GridView >(), std::declval< std::string >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        static const std::string clsName = name + std::to_string( dimRange );
        return FemPy::registerGridFunction< GridFunction >( scope, clsName.c_str() );
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridView, int dimRange >
      pybind11::object pyGlobalGridFunction ( const GridView &gridView, std::string name, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyGlobalGridFunction( gridView, std::move( name ), std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defGlobalGridFunction
    // ---------------------

    template< class GridView, int... dimRange >
    auto defGlobalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::registerPyGlobalGridFunction< GridView >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridView &, std::string, pybind11::function, pybind11::object ) > Dispatch;
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( detail::pyGlobalGridFunction< GridView, dimRange > )... }};

      return [ dispatch ] ( pybind11::object gridView, std::string name, pybind11::function evaluate ) {
          typename GridView::template Codim< 0 >::Geometry::GlobalCoordinate x( 0 );
          pybind11::gil_scoped_acquire acq;
          pybind11::object v( evaluate( x ) );
          const std::size_t dimR = len( v );
          if( dimR >= dispatch.size() )
            DUNE_THROW( NotImplemented, "globalGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
          return dispatch[ dimR ]( gridView.cast< const GridView & >(), std::move( name ), std::move( evaluate ), gridView );
        };
    }



    namespace detail
    {

      // makePyLocalGridFunction
      // -----------------------

      template< class GridView, int dimRange >
      auto makePyLocalGridFunction ( const GridView &gridView, std::string name, pybind11::function evaluate, std::integral_constant< int, dimRange > )
      {
        typedef typename GridView::template Codim< 0 >::Entity Entity;
        typedef typename GridView::template Codim< 0 >::Geometry::LocalCoordinate Coordinate;
        return simpleGridFunction( std::move( name ), gridView, [ evaluate ] ( const Entity &entity, const Coordinate &x ) {
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( entity, x ) );
            return v.template cast< FieldVector< double, dimRange > >();
          } );
      }


      // registerPyLocalGridFunction
      // ---------------------------

      template< class GridView, int dimRange >
      auto registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyLocalGridFunction( std::declval< GridView >(), std::declval< std::string >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        static const std::string clsName = name + std::to_string( dimRange );
        return FemPy::registerGridFunction< GridFunction >( scope, clsName.c_str() );
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridView, int dimRange >
      pybind11::object pyLocalGridFunction ( const GridView &gridView, std::string name, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyLocalGridFunction( gridView, std::move( name ), std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defLocalGridFunction
    // --------------------

    template< class GridView, int... dimRange >
    auto defLocalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::registerPyLocalGridFunction< GridView >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridView &, std::string, pybind11::function, pybind11::object ) > Dispatch;
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( detail::pyLocalGridFunction< GridView, dimRange > )... }};

      return [ dispatch ] ( pybind11::object gp, std::string name, pybind11::function evaluate ) {
          const GridView &gridView = gp.cast< const GridView & >();
          int dimR = -1;
          if( gridView.template begin< 0 >() != gridView.template end< 0 >() )
          {
            typename GridView::template Codim< 0 >::Geometry::LocalCoordinate x( 0 );
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( *gridView.template begin< 0 >(), x ) );
            dimR = len( v );
          }
          dimR = gridView.comm().max( dimR );
          if( dimR < 0 )
            DUNE_THROW( InvalidStateException, "Cannot create local grid function on empty grid" );
          if( static_cast< std::size_t >( dimR ) >= dispatch.size() )
            DUNE_THROW( NotImplemented, "localGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
          return dispatch[ static_cast< std::size_t >( dimR ) ]( gridView, std::move( name ), std::move( evaluate ), std::move( gp ) );
        };
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
