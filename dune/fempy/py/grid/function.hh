#ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
#define DUNE_FEMPY_PY_GRID_FUNCTION_HH

#include <functional>
#include <string>
#include <tuple>
#include <utility>

#include <dune/common/visibility.hh>

#if DUNE_VERSION_NEWER( DUNE_FEM, 2, 5 )
#include <dune/fem/misc/domainintegral.hh>
#else // #if DUNE_VERSION_NEWER( DUNE_FEM, 2, 5 )
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/grid/common/rangegenerators.hh>
namespace Dune {
  namespace Fem {
    template <class GridPart>
    struct Integral {
      typedef typename GridPart::template Codim< 0 >::EntityType ElementType;
      Integral(const GridPart &gp, int order)
        : gp_(gp), order_(order) {}
      template <class GF>
      typename GF::RangeType norm(const GF &gf) {
      typename GF::RangeType ret(0);
      const auto end = gp_.template end<0>();
      for (auto it = gp_.template begin<0>(); it!=end; ++it) {
          const auto &element = *it;
          auto lf = gf.localFunction(element);
          Dune::Fem::CachingQuadrature< GridPart, 0 > quadrature(element, order_);
          const size_t numQuadraturePoints = quadrature.nop();
          for( size_t pt = 0; pt < numQuadraturePoints; ++pt ) {
            const auto &x = quadrature.point( pt );
            const double weight = quadrature.weight( pt ) * element.geometry().integrationElement( x );
            typename GF::RangeType val;
            lf.evaluate( quadrature[ pt ], val );
            ret.axpy(weight,val);
          }
        }
        return ret;
      }
      private:
      const GridPart &gp_;
      int order_;
    };
  } // namespcae Fem
} // namespcae Dune
#endif // #else // #if DUNE_VERSION_NEWER( DUNE_FEM, 2, 5 )

#include <dune/corepy/grid/vtk.hh>

#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/numpy.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    // registerLocalFunction
    // ---------------------

    template< class LocalFunction >
    pybind11::class_< LocalFunction > registerLocalFunction ( pybind11::handle scope, const char *clsName = "LocalFunction" )
    {
      typedef typename LocalFunction::EntityType::Geometry::LocalCoordinate LocalCoordinate;

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

      // makePyLocalGridFunction
      // -----------------------

      template< class GridPart, int dimRange >
      auto makePyLocalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, std::integral_constant< int, dimRange > )
      {
        typedef typename GridPart::template Codim< 0 >::EntityType Entity;
        typedef typename GridPart::template Codim< 0 >::GeometryType::LocalCoordinate Coordinate;
        return simpleGridFunction( std::move( name ), gridPart, [ evaluate ] ( const Entity &entity, const Coordinate &x ) {
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( entity, x ) );
            try { return v.template cast< FieldVector< double, dimRange > >(); }
            catch(std::exception& e) { std::cout << e.what() << " in LocalGridFunction::evaluate" << std::endl; throw pybind11::cast_error("error converting return value in localGridFunction"); }
            return FieldVector<double,dimRange>(0);
          }, order );
      }

    }

    namespace detail
    {

      // registerGridFunction
      // --------------------
      template< class GridPart, int dimRange >
      auto registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > );

      template< class GridFunction, class Cls >
      void registerGridFunction ( pybind11::handle scope, Cls &cls)
      {
        using pybind11::operator""_a;

        typedef typename GridFunction::LocalFunctionType LocalFunction;
        typedef typename LocalFunction::EntityType Entity;
        typedef typename GridFunction::GridPartType GridPartType;
        typedef typename GridPartType::GridViewType GridView;

        registerLocalFunction< LocalFunction >( cls );

        cls.def( "__repr__", [] ( GridFunction &gf ) -> std::string {
            return "GridFunction< " + std::to_string( GridFunction::RangeType::dimension ) + " >(name = " + gf.name() + ")";
          } );

        cls.def_property_readonly( "dimRange", [] ( GridFunction &gf ) -> int { return GridFunction::RangeType::dimension; } );
        cls.def_property_readonly( "order", [] ( GridFunction &gf ) -> unsigned int { return gf.space().order(); } );
        cls.def_property_readonly( "name", [] ( GridFunction &gf ) -> std::string { return gf.name(); } );
        cls.def_property_readonly( "grid", [] ( GridFunction &gf ) -> GridView { return static_cast< GridView >( gf.gridPart() ); } );

        cls.def( "localFunction", [] ( const GridFunction &gf, const Entity &entity ) -> LocalFunction {
            return gf.localFunction( entity );
          }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >() );

        typedef decltype( makePyLocalGridFunction( std::declval< GridPartType >(), std::declval< std::string >(), std::declval< int >(), std::declval< pybind11::function >(), std::integral_constant< int, 1 >() ) ) ScalarLocalGridFunction;
        if (!pybind11::already_registered<ScalarLocalGridFunction>())
          registerPyLocalGridFunction<GridPartType>( scope, "LocalGridFunction", std::integral_constant< int, 1 >() );

        cls.def( "__getitem__", [] (const GridFunction &gf, size_t c ) {
            return makePyLocalGridFunction( gf.gridPart(), gf.name()+"_"+std::to_string(c), gf.space().order(),
                pybind11::cpp_function(
                [gf,c](const Entity &e, const typename Entity::Geometry::LocalCoordinate &x)
                {
                  typename GridFunction::LocalFunctionType::RangeType value;
                  gf.localFunction(e).evaluate( x, value );
                  return value[c];
                }), std::integral_constant< int, 1 >() );
          }, pybind11::keep_alive< 0, 1 >() );

        cls.def( "addToVTKWriter", &Dune::CorePy::addToVTKWriter< GridFunction >, pybind11::keep_alive< 3, 1 >(), "name"_a, "writer"_a, "dataType"_a );

        cls.def( "cellData", [] ( const GridFunction &self, int level ) { return cellData( self, level ); }, "level"_a = 0 );
        cls.def( "pointData", [] ( const GridFunction &self, int level ) { return pointData( self, level ); }, "level"_a = 0 );

        cls.def( "integrate", [] ( const GridFunction &gf ) { return Dune::Fem::Integral<GridPartType>(gf.gridPart(),gf.space().order()).norm(gf); });


#if 0
        cls.def_property_readonly( "as_ufl", [] ( GridFunction &gf ) -> pybind11::handle
            { pybind11::tuple args( 1 );
              args[ 0 ] = gf;
              return PyObject_Call( Dune::FemPy::getGridFunctionWrapper().ptr(), args.ptr(), nullptr );
            } );
#elif 0
        cls.def( "as_ufl", [] ( GridFunction &gf ) -> pybind11::handle
            {
              pybind11::tuple args( 1 );
              args[ 0 ] = gf;
              return PyObject_Call( Dune::FemPy::getGridFunctionWrapper().ptr(), args.ptr(), nullptr );
            }, pybind11::keep_alive<0,1>() );
#else
        cls.def( "as_ufl", [] ( pybind11::object &gf ) -> pybind11::handle
            { pybind11::tuple args( 1 );
              args[ 0 ] = gf;
              return PyObject_Call( Dune::FemPy::getGridFunctionWrapper().ptr(), args.ptr(), nullptr );
            },  pybind11::keep_alive<0,1>() );
#endif

      }

      template< class GridFunction >
      pybind11::class_< GridFunction > registerGridFunction ( pybind11::handle scope, const char *clsName = "GridFunction" )
      {
        pybind11::class_< GridFunction > cls( scope, clsName );
        registerGridFunction< GridFunction >( scope, cls );
        return cls;
      }


      // clsVirtualizedGridFunction
      // --------------------------

      template< class GridPart, class Value >
      DUNE_EXPORT inline pybind11::class_< VirtualizedGridFunction< GridPart, Value > > clsVirtualizedGridFunction ( pybind11::handle scope )
      {
        typedef VirtualizedGridFunction< GridPart, Value > GridFunction;
        static const std::string clsName = "VirtualizedGridFunction" + std::to_string( Value::dimension );
        static pybind11::class_< GridFunction > cls = registerGridFunction< GridFunction >( scope, clsName.c_str() );
        return cls;
      }

    } // namespace detail



    // registerGridFunction
    // --------------------

    template< class GridFunction >
    pybind11::class_< GridFunction > registerGridFunction ( pybind11::handle scope, const char *clsName = "GridFunction" )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Value;

      auto cls = detail::registerGridFunction< GridFunction >( scope, clsName );

      detail::clsVirtualizedGridFunction< GridPart, Value >( scope ).def( pybind11::init< GridFunction >() );
      pybind11::implicitly_convertible< GridFunction, VirtualizedGridFunction< GridPart, Value > >();

      return cls;
    }



    // registerVirtualizedGridFunction
    // -------------------------------

    template< class GridPart, int... dimRange >
    void registerVirtualizedGridFunction ( pybind11::handle scope, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::clsVirtualizedGridFunction< GridPart, FieldVector< double, dimRange > >( scope )... );
    };



    namespace detail
    {

      // makePyGlobalGridFunction
      // ------------------------

      template< class GridPart, int dimRange >
      auto makePyGlobalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, std::integral_constant< int, dimRange > )
      {
        typedef typename GridPart::template Codim< 0 >::GeometryType::GlobalCoordinate Coordinate;
        return simpleGridFunction( std::move( name ), gridPart, [ evaluate ] ( const Coordinate &x ) {
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( x ) );
            try { return v.template cast< FieldVector< double, dimRange > >(); }
            catch(std::exception& e) { std::cout << e.what() << " in GlobalGridFunction::evaluate" << std::endl; throw pybind11::cast_error("error converting return value in localGridFunction"); }
            return FieldVector<double,dimRange>(0);
          }, order );
      }



      // registerPyGlobalGridFunction
      // ----------------------------

      template< class GridPart, int dimRange >
      auto registerPyGlobalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyGlobalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< int >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        static const std::string clsName = name + std::to_string( dimRange );
        return FemPy::registerGridFunction< GridFunction >( scope, clsName.c_str() );
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridPart, int dimRange >
      pybind11::object pyGlobalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyGlobalGridFunction( gridPart, std::move( name ), order, std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defGlobalGridFunction
    // ---------------------

    template< class GridPart, int... dimRange >
    auto defGlobalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::registerPyGlobalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridPart &, std::string, int, pybind11::function, pybind11::object ) > Dispatch;
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( detail::pyGlobalGridFunction< GridPart, dimRange > )... }};

      return [ dispatch ] ( pybind11::object gv, std::string name, int order, pybind11::function evaluate ) {
          typedef typename GridPart::GridViewType GridView;
          const auto &gp = gridPart<GridView>( gv );
          // typename GridPart::template Codim< 0 >::GeometryType::GlobalCoordinate x( 0 );
          auto x = gp.template begin<0>()->geometry().center();
          pybind11::gil_scoped_acquire acq;
          pybind11::object v( evaluate( x ) );
          const std::size_t dimR = len( v );
          if( dimR >= dispatch.size() )
            DUNE_THROW( NotImplemented, "globalGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
          return dispatch[ dimR ]( gp, std::move( name ), order, std::move( evaluate ), std::move(gv) );
        };
    }



    namespace detail
    {

      // registerPyLocalGridFunction
      // ---------------------------

      template< class GridPart, int dimRange >
      DUNE_EXPORT auto registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyLocalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< int >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        static const std::string clsName = name + std::to_string( dimRange );
        static const auto cls = FemPy::registerGridFunction< GridFunction >( scope, clsName.c_str() );
        return cls;
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridPart, int dimRange >
      pybind11::object pyLocalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyLocalGridFunction( gridPart, std::move( name ), order, std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defLocalGridFunction
    // --------------------

    template< class GridPart, int... dimRange >
    auto defLocalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::registerPyLocalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridPart &, std::string, int, pybind11::function, pybind11::object ) > Dispatch;
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( detail::pyLocalGridFunction< GridPart, dimRange > )... }};

      return [ dispatch ] ( pybind11::object gv, std::string name, int order, pybind11::function evaluate ) {
          typedef typename GridPart::GridViewType GridView;
          const auto &gp = gridPart<GridView>( gv );
          int dimR = -1;
          if( gp.template begin< 0 >() != gp.template end< 0 >() )
          {
            typename GridPart::template Codim< 0 >::GeometryType::LocalCoordinate x( 0 );
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( *gp.template begin< 0 >(), x ) );
            dimR = len( v );
          }
          dimR = gp.comm().max( dimR );
          if( dimR < 0 )
            DUNE_THROW( InvalidStateException, "Cannot create local grid function on empty grid" );
          if( static_cast< std::size_t >( dimR ) >= dispatch.size() )
            DUNE_THROW( NotImplemented, "localGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
          return dispatch[ static_cast< std::size_t >( dimR ) ]( gp, std::move( name ), order, std::move( evaluate ), std::move( gv ) );
        };
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
