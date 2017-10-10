#ifndef DUNE_FEMPY_FUNCTION_SUBGRIDFUNCTION_HH
#define DUNE_FEMPY_FUNCTION_SUBGRIDFUNCTION_HH

#include <cstddef>

#include <functional>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/fempy/function/utility.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      template< class FunctionSpace >
      struct ScalarFunctionSpace;

      template< class DF, class RF, int dimD, int dimR >
      struct ScalarFunctionSpace< Fem::FunctionSpace< DF, RF, dimD, dimR > >
      {
        typedef Fem::FunctionSpace< DF, RF, dimD, 1 > Type;
      };

    } // namespace detail


    template< class FunctionSpace >
    using ScalarFunctionSpaceType = typename detail::ScalarFunctionSpace< FunctionSpace >::Type;



    // SubLocalFunction
    // ----------------

    template< class LocalFunction >
    class SubLocalFunction
    {
      typedef SubLocalFunction< LocalFunction > This;

    public:
      typedef typename LocalFunction::EntityType EntityType;

      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;
      typedef typename EntityType::Geometry::GlobalCoordinate GlobalCoordinateType;

      typedef ScalarFunctionSpaceType< typename LocalFunction::FunctionSpaceType > FunctionSpaceType;

      static const int dimDomain = FunctionSpaceType::dimDomain;
      static const int dimRange = FunctionSpaceType::dimRange;

      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      SubLocalFunction ( LocalFunction &&localFunction, std::size_t component )
        : localFunction_( std::move( localFunction ) ), component_( component )
      {}

      template< class GF,
                std::enable_if_t< std::is_constructible< LocalFunction, decltype( std::declval< const GF & >().gridFunction() ) >::value, int > = 0,
                std::enable_if_t< std::is_constructible< std::size_t, decltype( std::declval< const GF & >().component() ) >::value, int > = 0 >
      explicit SubLocalFunction ( const GF &gf )
        : localFunction_( gf.gridFunction() ), component_( gf.component() )
      {}

      template< class GF,
                std::enable_if_t< std::is_constructible< LocalFunction, decltype( std::declval< const GF & >().gridFunction() ) >::value, int > = 0,
                std::enable_if_t< std::is_constructible< std::size_t, decltype( std::declval< const GF & >().component() ) >::value, int > = 0 >
      explicit SubLocalFunction ( const GF &gf, const EntityType &entity )
        : localFunction_( gf.gridFunction() ), component_( gf.component() )
      {
        init( entity );
      }

      void init ( const EntityType &entity ) { localFunction().init( entity ); }

      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        typename LocalFunction::RangeType tmp;
        localFunction().evaluate( x, tmp );
        value[ 0 ] = tmp[ component() ];
      }

      template< class Quadrature, class Values >
      void evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
      {
        for( const auto qp : quadrature )
          evaluate( qp, values[ qp.index() ] );
      }

      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        typename LocalFunction::JacobianRangeType tmp;
        localFunction().jacobian( x, tmp );
        jacobian[ 0 ] = tmp[ component() ];
      }

      template< class Quadrature, class Jacobians >
      void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
      {
        for( const auto qp : quadrature )
          jacobian( qp, jacobians[ qp.index() ] );
      }

      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        typename LocalFunction::HessianRangeType tmp;
        localFunction().hessian( x, tmp );
        hessian[ 0 ] = tmp[ component() ];
      }

      int order () const { return localFunction().order(); }

      const EntityType &entity () const { return localFunction().entity(); }

      const LocalFunction &localFunction () const { return localFunction_; }
      LocalFunction &localFunction () { return localFunction_; }

      std::size_t component () const { return component_; }

    private:
      LocalFunction localFunction_;
      std::size_t component_;
    };



    // SubGridFunction
    // ---------------

    template< class GridFunction >
    class SubGridFunction
      : public Fem::Function< ScalarFunctionSpaceType< typename std::decay_t< decltype( std::ref( std::declval< GridFunction & >() ).get() ) >::FunctionSpaceType >, SubGridFunction< GridFunction > >,
        public Fem::HasLocalFunction
    {
      typedef SubGridFunction< GridFunction > This;
      typedef Fem::Function< ScalarFunctionSpaceType< typename std::decay_t< decltype( std::ref( std::declval< GridFunction & >() ).get() ) >::FunctionSpaceType >, SubGridFunction< GridFunction > > Base;

      static_assert( std::is_same< ScalarFunctionSpaceType< typename Base::FunctionSpaceType >, typename Base::FunctionSpaceType >::value, "SubGridFunction does not have a scalar function space." );

    public:
      typedef std::decay_t< decltype( std::ref( std::declval< GridFunction & >() ).get() ) > GridFunctionType;

      typedef typename GridFunctionType::GridPartType GridPartType;

      typedef SubLocalFunction< Fem::ConstLocalFunction< GridFunctionType > > LocalFunctionType;

      typedef typename LocalFunctionType::EntityType EntityType;

      using typename Base::FunctionSpaceType;

      static const int dimDomain = FunctionSpaceType::dimDomain;
      static const int dimRange = FunctionSpaceType::dimRange;

      using typename Base::DomainType;
      using typename Base::RangeType;
      using typename Base::JacobianRangeType;

      SubGridFunction ( GridFunction gridFunction, std::size_t component )
        : gridFunction_( std::move( gridFunction ) ), component_( component )
      {}

      LocalFunctionType localFunction () const { return LocalFunctionType( *this ); }
      LocalFunctionType localFunction ( const EntityType &entity ) const { return LocalFunctionType( *this, entity ); }

      std::string name () const { return gridFunction().name() + "[ " + std::to_string( component() ) + " ]"; }

      const GridPartType &gridPart () const { return gridFunction().gridPart(); }

      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        typename GridFunction::RangeType tmp;
        gridFunction().evaluate( x, tmp );
        value[ 0 ] = tmp[ component_ ];
      }

      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        typename GridFunction::JacobianRangeType tmp;
        gridFunction().evaluate( x, tmp );
        jacobian[ 0 ] = jacobian[ tmp ];
      }

      int order () const { return FemPy::order( gridFunction() ); }

      const GridFunctionType &gridFunction () const { return std::ref( gridFunction_ ).get(); }
      const std::size_t component () const { return component_; }

    protected:
      GridFunction gridFunction_;
      std::size_t component_;
    };



    // subGridFunction
    // ---------------

    template< class GridFunction, std::enable_if_t< !std::is_lvalue_reference< GridFunction >::value, int > = 0 >
    inline static SubGridFunction< GridFunction > subGridFunction ( GridFunction &&gridFunction, std::size_t component )
    {
      return SubGridFunction< GridFunction >( std::move( gridFunction ), component );
    }

    template< class GridFunction, std::enable_if_t< std::is_lvalue_reference< GridFunction >::value, int > = 0 >
    inline static SubGridFunction< std::reference_wrapper< const GridFunction > > subGridFunction ( GridFunction &&gridFunction, const std::size_t component )
    {
      return SubGridFunction< std::reference_wrapper< const GridFunction > >( std::ref( gridFunction ), component );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_SUBGRIDFUNCTION_HH
