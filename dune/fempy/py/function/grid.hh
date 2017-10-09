#ifndef DUNE_FEMPY_PY_FUNCTION_GRID_HH
#define DUNE_FEMPY_PY_FUNCTION_GRID_HH

#include <functional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/classname.hh>
#include <dune/common/visibility.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/misc/domainintegral.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/corepy/common/string_constant.hh>
#include <dune/corepy/grid/vtk.hh>

#include <dune/fempy/common/utility.hh>
#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fempy/function/subgridfunction.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/numpy.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class GridFunction, class... options >
    static void registerGridFunction ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls );



    // PyGridFunction
    // --------------

    template< class GridFunction >
    class DUNE_PRIVATE PyGridFunction
      : public Fem::Function< typename GridFunction::FunctionSpaceType, PyGridFunction< GridFunction > >,
        public Fem::HasLocalFunction
    {
      typedef Fem::Function< typename GridFunction::FunctionSpaceType, PyGridFunction< GridFunction > > Base;

    public:
      typedef typename GridFunction::GridPartType GridPartType;

      typedef typename GridFunction::EntityType EntityType;

      typedef typename Base::DomainType DomainType;
      typedef typename Base::RangeType RangeType;
      typedef typename Base::JacobianRangeType JacobianRangeType;
      typedef typename Base::HessianRangeType HessianRangeType;

      class LocalFunctionType
      {
        typedef typename GridFunction::LocalFunctionType Impl;

      public:
        typedef typename Impl::EntityType EntityType;

        typedef typename Impl::FunctionSpaceType FunctionSpaceType;

        static const int dimDomain = FunctionSpaceType::dimDomain;
        static const int dimRange = FunctionSpaceType::dimRange;

        typedef typename FunctionSpaceType::DomainType DomainType;
        typedef typename FunctionSpaceType::RangeType RangeType;
        typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

        explicit LocalFunctionType ( const PyGridFunction &gf ) : impl_( *gf.impl_ ), pyObj_( gf.pyObj_ ) {}

        LocalFunctionType ( Impl impl, pybind11::object pyObj ) : impl_( std::move( impl ) ), pyObj_( std::move( pyObj ) ) {}

        void init ( const EntityType &entity ) { impl_.init( entity ); }

        template< class Point >
        void evaluate ( const Point &x, RangeType &value ) const
        {
          impl_.evaluate( x, value );
        }

        template< class Quadrature, class Values >
        void evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
        {
          impl_.evaluateQuadrature( quadrature, values );
        }

        template< class Point >
        void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
        {
          impl_.jacobian( x, jacobian );
        }

        template< class Quadrature, class Jacobians >
        void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
        {
          impl_.jacobianQuadrature( quadrature, jacobians );
        }

        template< class Point >
        void hessian ( const Point &x, HessianRangeType &hessian ) const
        {
          impl_.hessian( x, hessian );
        }

        int order () const { return impl_.order(); }

        const EntityType &entity () const { return impl_.entity(); }

      private:
        Impl impl_;
        pybind11::object pyObj_;
      };

    public:
      PyGridFunction ( const GridFunction &impl, pybind11::object pyObj )
        : impl_( &impl ), pyObj_( std::move( pyObj ) )
      {}

      PyGridFunction ( const GridFunction &impl )
        : impl_( &impl ),
          pyObj_( pybind11::detail::get_object_handle( impl_, pybind11::detail::get_type_info( typeid( GridFunction ) ) ), true )
      {}

      LocalFunctionType localFunction ( const EntityType &entity ) const { return LocalFunctionType( impl_->localFunction( entity ), pyObj_ ); }

      std::string name () const { return impl_->name(); }

      const GridPartType &gridPart () const { return impl_->gridPart(); }

      // !!!! void evaluate ( const DomainType &x, RangeType &value ) const { return impl_->evaluate( x, value ); }
      // !!!! void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const { return impl_->jacobian( x, jacobian ); }

    protected:
      const GridFunction *impl_;
      pybind11::object pyObj_;
    };



    // pyGridFunction
    // --------------

    template< class GridFunction >
    inline static PyGridFunction< GridFunction > pyGridFunction ( const GridFunction &gridFunction ) noexcept
    {
      return PyGridFunction< GridFunction >( gridFunction );
    }

    template< class GridFunction >
    inline static PyGridFunction< GridFunction > pyGridFunction ( const GridFunction &gridFunction, pybind11::object pyObj ) noexcept
    {
      return PyGridFunction< GridFunction >( gridFunction, std::move( pyObj ) );
    }

    template< class GridFunction >
    inline static PyGridFunction< GridFunction > pyGridFunction ( pybind11::object pyObj ) noexcept
    {
      return PyGridFunction< GridFunction >( pybind11::cast< const GridFunction & >( pyObj ), std::move( pyObj ) );
    }



    // registerLocalFunction
    // ---------------------

    template< class LocalFunction, class... options >
    inline static void registerLocalFunction ( pybind11::handle scope, pybind11::class_< LocalFunction, options... > cls )
    {
      typedef typename LocalFunction::EntityType::Geometry::LocalCoordinate LocalCoordinate;

      using pybind11::operator""_a;

      cls.def_property_readonly( "dimRange", [] ( LocalFunction & ) -> int { return LocalFunction::RangeType::dimension; } );
      cls.def( "evaluate", [] ( const LocalFunction &self, const LocalCoordinate &x ) {
          typename LocalFunction::RangeType value;
          self.evaluate( x, value );
          return value;
        }, "x"_a );
      cls.def( "jacobian", [] ( const LocalFunction &self, const LocalCoordinate &x ) {
          typename LocalFunction::JacobianRangeType jacobian;
          self.jacobian( x, jacobian );
          return jacobian;
        }, "x"_a );
    }

    template< class LocalFunction >
    inline static pybind11::class_< LocalFunction > registerLocalFunction ( pybind11::handle scope, const char *clsName = "LocalFunction" )
    {
      pybind11::class_< LocalFunction > cls( scope, clsName );
      registerLocalFunction( scope, cls );
      return cls;
    }



    namespace detail
    {

      // registerGridFunctionSubscript
      // -----------------------------

      template< class GridFunction, class... options >
      inline static auto registerGridFunctionSubscript ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< !std::is_same< ScalarFunctionSpace< typename GridFunction::FunctionSpaceType >, typename GridFunction::FunctionSpaceType >::value >
      {
        typedef FemPy::SubGridFunction< PyGridFunction< GridFunction > > SubGridFunction;

        using pybind11::operator""_a;

        FemPy::registerGridFunction< SubGridFunction >( cls, "SubFunction" );

        cls.def( "__getitem__", [] ( pybind11::object self, std::size_t c ) {
            if( c < static_cast< std::size_t >( GridFunction::RangeType::dimension ) )
              return SubGridFunction( pyGridFunction< GridFunction >( self ), c );
            else
              throw pybind11::index_error();
          }, "c"_a );
      }

      template< class GridFunction, class... options >
      inline static auto registerGridFunctionSubscript ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_same< ScalarFunctionSpace< typename GridFunction::FunctionSpaceType >, typename GridFunction::FunctionSpaceType >::value >
      {
        using pybind11::operator""_a;

        cls.def( "__getitem__", [] ( pybind11::object self, std::size_t c ) {
            if( c == 0 )
              return self;
            else
              throw pybind11::index_error();
          }, "c"_a );
      }

      template< class GridFunction, class... options >
      inline static void registerGridFunctionSubscript ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls, PriorityTag< 0 > )
      {}

      template< class GridFunction, class... options >
      inline static void registerGridFunctionSubscript ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls )
      {
        registerGridFunctionSubscript( scope, cls, PriorityTag< 42 >() );
      }



      // registerGridFunction
      // --------------------

      template< class GridFunction, class... options >
      inline static void registerGridFunction ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls )
      {
        using pybind11::operator""_a;

        typedef typename GridFunction::LocalFunctionType LocalFunction;
        typedef typename LocalFunction::EntityType Entity;
        typedef typename GridFunction::GridPartType GridPartType;
        typedef typename GridPartType::GridViewType GridView;

        registerLocalFunction< LocalFunction >( cls );

        cls.def( "__repr__", [] ( GridFunction &self ) -> std::string {
            return "GridFunction< " + std::to_string( GridFunction::RangeType::dimension ) + " >(name = " + self.name() + ")";
          } );

        cls.def_property_readonly( "dimRange", [] ( GridFunction &self ) -> int { return GridFunction::RangeType::dimension; } );
        cls.def_property_readonly( "order", [] ( GridFunction &self ) -> unsigned int { return self.space().order(); } );
        cls.def_property_readonly( "name", [] ( GridFunction &self ) -> std::string { return self.name(); } );
        cls.def_property_readonly( "grid", [] ( GridFunction &self ) -> GridView { return static_cast< GridView >( self.gridPart() ); } );

        cls.def( "localFunction", [] ( const GridFunction &self, const Entity &entity ) -> LocalFunction {
            return self.localFunction( entity );
          }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >() );

        cls.def( "addToVTKWriter", &Dune::CorePy::addToVTKWriter< GridFunction >, pybind11::keep_alive< 3, 1 >(), "name"_a, "writer"_a, "dataType"_a );

        cls.def( "cellData", [] ( const GridFunction &self, int level ) { return cellData( self, refinementLevels( level ) ); }, "level"_a = 0 );
        cls.def( "pointData", [] ( const GridFunction &self, int level ) { return pointData( self, refinementLevels( level ) ); }, "level"_a = 0 );

        cls.def( "integrate", [] ( const GridFunction &self ) { return Dune::Fem::Integral<GridPartType>(self.gridPart(),self.space().order()).norm(self); });
      }



      // nameVirtualizedGridFunction
      // ---------------------------

      template< int dimRange >
      static const auto nameVirtualizedGridFunction = "VirtualizedGridFunction_" + CorePy::make_string_constant( std::integral_constant< int, dimRange >() );



      // registerVirtualizedGridFunction
      // -------------------------------

      template< class GridPart, class Value >
      inline void registerVirtualizedGridFunction ( pybind11::handle scope )
      {
        typedef VirtualizedGridFunction< GridPart, Value > GridFunction;
        if( !pybind11::already_registered< GridFunction >() )
          detail::registerGridFunction( scope, pybind11::class_< GridFunction >( scope, nameVirtualizedGridFunction< Value::dimension > ) );
      }

    } // namespace detail



    // registerGridFunction
    // --------------------

    template< class GridFunction, class... options >
    inline static void registerGridFunction ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Value;

      detail::registerGridFunction( scope, cls );

      detail::registerVirtualizedGridFunction< GridPart, Value >( scope );
      cls.def( "asVirtualizedGridFunction", [] ( GridFunction &self ) { return new VirtualizedGridFunction< GridPart, Value >( pyGridFunction( self ) ); } );
    }

    template< class GridFunction >
    inline static pybind11::class_< GridFunction > registerGridFunction ( pybind11::handle scope, const char *clsName = "GridFunction" )
    {
      pybind11::class_< GridFunction > cls( scope, clsName );
      FemPy::registerGridFunction( scope, cls );
      return cls;
    }



    namespace detail
    {

      template< class GF >
      using VirtualizedGridFunctionFor = VirtualizedGridFunction< typename GF::GridPartType, typename GF::RangeType >;



      // asGridFunction
      // --------------

      template< class GridFunction, class Apply >
      inline static auto asGridFunction ( const typename GridFunction::GridPartType &gridPart, pybind11::object gf, Apply &&apply, PriorityTag< 2 > )
        -> std::enable_if_t< std::is_same< GridFunction, VirtualizedGridFunctionFor< GridFunction > >::value, decltype( apply( std::declval< const GridFunction & >() ) ) >
      {
        typedef typename GridFunction::GridPartType::template Codim< 0 >::EntityType Entity;
        typedef typename GridFunction::GridPartType::template Codim< 0 >::GeometryType::LocalCoordinate LocalCoordinate;

        try
        {
          return apply( pybind11::cast< GridFunction >( gf ) );
        }
        catch( pybind11::cast_error )
        {}

        try
        {
          return apply( pybind11::cast< GridFunction >( gf.attr( "asVirtualizedGridFunction" )() ) );
        }
        catch( pybind11::error_already_set )
        {}

        try
        {
          const typename GridFunction::RangeType value = pybind11::cast< typename GridFunction::RangeType >( gf );
          return apply( simpleGridFunction( gridPart, [ value ] ( const Entity &entity, const LocalCoordinate & ) { return value; }, 0 ) );
        }
        catch( pybind11::cast_error )
        {}

        const std::string typeName( pybind11::str( gf.get_type() ) );
        throw pybind11::cast_error("Unable to cast Python instance of type " + typeName + " to grid function");
      }

      template< class GridFunction, class Apply >
      inline static auto asGridFunction ( const typename GridFunction::GridPartType &gridPart, pybind11::object gf, Apply &&apply, PriorityTag< 1 > )
        -> same_type< decltype( apply( std::declval< const GridFunction & >() ) ), decltype( apply( std::declval< const VirtualizedGridFunctionFor< GridFunction > & >() ) ) >
      {
        try
        {
          return apply( pybind11::cast< GridFunction >( gf ) );
        }
        catch( pybind11::cast_error )
        {}

        return detail::asGridFunction< VirtualizedGridFunctionFor< GridFunction > >( gridPart, gf, apply, PriorityTag< 42 >() );
      }

      template< class GridFunction, class Apply >
      inline static auto asGridFunction ( const typename GridFunction::GridPartType &gridPart, pybind11::object gf, Apply &&apply, PriorityTag< 0 > )
        -> decltype( apply( std::declval< const GridFunction & >() ) )
      {
        try
        {
          return apply( pybind11::cast< GridFunction >( gf ) );
        }
        catch( pybind11::cast_error )
        {}

        const std::string typeName( pybind11::str( gf.get_type() ) );
        throw pybind11::cast_error("Unable to cast Python instance of type " + typeName + " to grid function of fixed type '" + className< GridFunction >() + "'");
      }

    } // namespace detail



    // asGridFunction
    // --------------

    template< class GridFunction, class Apply >
    inline static auto asGridFunction ( const typename GridFunction::GridPartType &gridPart, pybind11::object gf, Apply &&apply )
      -> decltype( apply( std::declval< const GridFunction & >() ) )
    {
      return detail::asGridFunction< GridFunction >( gridPart, gf, apply, PriorityTag< 42 >() );
    }



    namespace detail
    {

      // makePyGlobalGridFunction
      // ------------------------

      template< class GridPart, int dimRange >
      inline static auto makePyGlobalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, std::integral_constant< int, dimRange > )
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
      inline static auto registerPyGlobalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyGlobalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< int >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        static const std::string clsName = name + std::to_string( dimRange );
        return FemPy::registerGridFunction< GridFunction >( scope, clsName.c_str() );
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridPart, int dimRange >
      inline static pybind11::object pyGlobalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyGlobalGridFunction( gridPart, std::move( name ), order, std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defGlobalGridFunction
    // ---------------------

    template< class GridPart, int... dimRange >
    inline static auto defGlobalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
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

      // makePyLocalGridFunction
      // -----------------------

      template< class GridPart, int dimRange >
      inline static auto makePyLocalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, std::integral_constant< int, dimRange > )
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



      // registerPyLocalGridFunction
      // ---------------------------

      template< class GridPart, int dimRange >
      DUNE_EXPORT inline auto registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyLocalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< int >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        static const std::string clsName = name + std::to_string( dimRange );
        static const auto cls = FemPy::registerGridFunction< GridFunction >( scope, clsName.c_str() );
        return cls;
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridPart, int dimRange >
      inline static pybind11::object pyLocalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyLocalGridFunction( gridPart, std::move( name ), order, std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defLocalGridFunction
    // --------------------

    template< class GridPart, int... dimRange >
    inline static auto defLocalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
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

#endif // #ifndef DUNE_FEMPY_PY_FUNCTION_GRID_HH
