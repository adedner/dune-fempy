#ifndef DUNE_FEMPY_FUNCTION_SIMPLEGRIDFUNCTION_HH
#define DUNE_FEMPY_FUNCTION_SIMPLEGRIDFUNCTION_HH

#include <cassert>

#include <limits>
#include <type_traits>
#include <utility>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>

namespace Dune
{

  namespace FemPy
  {

    //! \brief Traits class for vector function spaces
    template< class DomainField, class RangeField, int dimD, int dimR>
    struct FunctionSpace
    {
      /** \copydoc Dune::Fem::FunctionSpaceInterface::DomainFieldType */
      typedef DomainField DomainFieldType;
      /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeFieldType */
      typedef RangeField RangeFieldType;

      /** \brief dimension of range vector space */
      enum { dimDomain = dimD };
      /** \brief dimension of domain vector space */
      enum { dimRange = dimR };

      /** \copydoc Dune::Fem::FunctionSpaceInterface::DomainType */
      typedef FieldVector< DomainFieldType, dimDomain > DomainType;
      /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeType */
      typedef FieldVector< RangeFieldType, dimRange> RangeType;
      /** \brief linear mapping type */
      typedef FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType;
    };
    // SimpleLocalFunction
    // -------------------

    template< class GridView, class LocalEvaluator >
    class SimpleLocalFunction
    {
      typedef SimpleLocalFunction< GridView, LocalEvaluator > This;

    public:
      typedef typename GridView::template Codim< 0 >::Entity EntityType;

      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;
      typedef std::decay_t< std::result_of_t< LocalEvaluator( EntityType, LocalCoordinateType ) > > Value;

      typedef FunctionSpace< typename GridView::ctype, typename FieldTraits< Value >::field_type, GridView::dimensionworld, Value::dimension > FunctionSpaceType;

      static const int dimDomain = FunctionSpaceType::dimDomain;
      static const int dimRange = FunctionSpaceType::dimRange;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

      template< class GridFunction, std::enable_if_t< std::is_same< This, typename GridFunction::LocalFunctionType >::value, int > = 0 >
      SimpleLocalFunction ( const GridFunction &gridFunction )
        : localEvaluator_( gridFunction.localEvaluator() ), order_( gridFunction.order() )
      {}

      SimpleLocalFunction ( const EntityType &entity, LocalEvaluator localEvaluator, int order )
        : entity_( &entity ), localEvaluator_( std::move( localEvaluator ) ), order_( order )
      {}

      void init ( const EntityType &entity ) { entity_ = &entity; }

      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        value = localEvaluator_( entity(), x );
      }

      template< class Quadrature, class Values >
      void evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
      {
        for( const auto qp : quadrature )
          evaluate( qp, values[ qp.index() ] );
      }

      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        DUNE_THROW( NotImplemented, "SimpleLocalFunction::jacobian not implemented" );
      }

      template< class Quadrature, class Jacobians >
      void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
      {
        for( const auto qp : quadrature )
          jacobian( qp, jacobians[ qp.index() ] );
      }

      int order () const { return order_; }

      const EntityType &entity () const { assert( entity_ ); return *entity_; }

    private:
      const EntityType *entity_ = nullptr;
      LocalEvaluator localEvaluator_;
      int order_;
    };



    // IsLocalEvaluator
    // ----------------

    template< class LE, class E, class X >
    std::true_type __isLocalEvaluator ( const LE &, const E &, const X &, decltype( std::declval< LE >()( std::declval< E >(), std::declval< X >() ) ) * = nullptr );

    std::false_type __isLocalEvaluator ( ... );

    template< class GP, class LE >
    struct IsLocalEvaluator
      : public decltype( __isLocalEvaluator( std::declval< LE >(), std::declval< typename GP::template Codim< 0 >::Entity >(), std::declval< typename GP::template Codim< 0 >::Geometry::LocalCoordinate >() ) )
    {};



    // SimpleGridFunction
    // ------------------

    template< class GridView, class LocalEvaluator >
    class SimpleGridFunction
    {
      typedef SimpleGridFunction< GridView, LocalEvaluator > This;

    public:
      typedef GridView GridViewType;

      typedef SimpleLocalFunction< GridView, LocalEvaluator > LocalFunctionType;

      typedef typename LocalFunctionType::EntityType EntityType;

      typedef typename LocalFunctionType::FunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

      SimpleGridFunction ( std::string name, const GridViewType &gridView, LocalEvaluator localEvaluator, int order = std::numeric_limits< int >::max() )
        : name_( std::move( name ) ),
          gridView_( gridView ),
          localEvaluator_( std::move( localEvaluator ) ),
          order_( order )
      {}

      SimpleGridFunction ( const GridViewType &gridView, LocalEvaluator localEvaluator, int order = std::numeric_limits< int >::max() )
        : name_( "unnamed-" + className< LocalEvaluator >() ),
          gridView_( gridView ),
          localEvaluator_( std::move( localEvaluator ) ),
          order_( order )
      {}

      LocalFunctionType localFunction ( const EntityType &entity ) const { return LocalFunctionType( entity, localEvaluator_, order_ ); }

      const std::string &name () const { return name_; }

      const GridViewType &gridView () const { return gridView_; }

      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        // needs to be implemented
      }

      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        // needs to be implemented
      }

      int order () const { return order_; }

      const LocalEvaluator &localEvaluator () const { return localEvaluator_; }

    protected:
      std::string name_;
      const GridViewType &gridView_;
      LocalEvaluator localEvaluator_;
      int order_;
    };



    // IsEvaluator
    // -----------

    template< class E, class X >
    std::true_type __isEvaluator ( const E &, const X &, decltype( std::declval< E >()( std::declval< X >() ) ) * = nullptr );

    std::false_type __isEvaluator ( ... );

    template< class GP, class E >
    struct IsEvaluator
      : public decltype( __isEvaluator( std::declval< E >(), std::declval< typename GP::template Codim< 0 >::Geometry::GlobalCoordinate >() ) )
    {};



    // LocalEvaluatorAdapter
    // ---------------------

    template< class Entity, class Evaluator >
    struct LocalEvaluatorAdapter
    {
      typedef typename Entity::Geometry::GlobalCoordinate GlobalCoordinate;
      typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;

      typedef decltype( std::declval< Evaluator >()( std::declval< GlobalCoordinate >() ) ) Value;

      LocalEvaluatorAdapter ( Evaluator evaluator ) : evaluator_( std::move( evaluator ) ) {}

      Value operator () ( const GlobalCoordinate &x ) const { return evaluator_( x ); }
      Value operator () ( const Entity &entity, const LocalCoordinate &x ) const { return evaluator_( entity.geometry().global( x ) ); }

    private:
      Evaluator evaluator_;
    };



    // SimpleGlobalGridFunction
    // ------------------------

    template< class GridView, class Evaluator >
    class SimpleGlobalGridFunction
      : public SimpleGridFunction< GridView, LocalEvaluatorAdapter< typename GridView::template Codim< 0 >::Entity, Evaluator > >
    {
      typedef SimpleGlobalGridFunction< GridView, Evaluator > This;
      typedef SimpleGridFunction< GridView, LocalEvaluatorAdapter< typename GridView::template Codim< 0 >::Entity, Evaluator > > Base;

    public:
      typedef typename Base::GridViewType GridViewType;

      typedef typename Base::DomainType DomainType;
      typedef typename Base::RangeType RangeType;

      SimpleGlobalGridFunction ( std::string name, const GridViewType &gridView, Evaluator evaluator, int order = std::numeric_limits< int >::max() )
        : Base( std::move( name ), gridView, std::move( evaluator ), order )
      {}

      SimpleGlobalGridFunction ( const GridViewType &gridView, Evaluator evaluator, int order = std::numeric_limits< int >::max() )
        : Base( "unnamed-" + className< Evaluator >(), gridView, std::move( evaluator ), order )
      {}

      void evaluate ( const DomainType &x, RangeType &value ) const { value = localEvaluator_( x ); }

    protected:
      using Base::localEvaluator_;
    };



    // simpleGridFunction
    // ------------------

    template< class GridView, class LocalEvaluator >
    inline static std::enable_if_t< IsLocalEvaluator< GridView, LocalEvaluator >::value, SimpleGridFunction< GridView, LocalEvaluator > >
    simpleGridFunction ( std::string name, const GridView &gridView, LocalEvaluator localEvaluator, int order = std::numeric_limits< int >::max() )
    {
      return SimpleGridFunction< GridView, LocalEvaluator >( std::move( name ), gridView, std::move( localEvaluator ), order );
    }

    template< class GridView, class LocalEvaluator >
    inline static std::enable_if_t< IsLocalEvaluator< GridView, LocalEvaluator >::value, SimpleGridFunction< GridView, LocalEvaluator > >
    simpleGridFunction ( const GridView &gridView, LocalEvaluator localEvaluator, int order = std::numeric_limits< int >::max() )
    {
      return SimpleGridFunction< GridView, LocalEvaluator >( gridView, std::move( localEvaluator ), order );
    }


    // simpleGridFunction
    // ------------------

    template< class GridView, class Evaluator >
    inline static std::enable_if_t< IsEvaluator< GridView, Evaluator >::value, SimpleGlobalGridFunction< GridView, Evaluator > >
    simpleGridFunction ( std::string name, const GridView &gridView, Evaluator evaluator, int order = std::numeric_limits< int >::max() )
    {
      return SimpleGlobalGridFunction< GridView, Evaluator >( std::move( name ), gridView, std::move( evaluator ), order );
    }

    template< class GridView, class Evaluator >
    inline static std::enable_if_t< IsEvaluator< GridView, Evaluator >::value, SimpleGlobalGridFunction< GridView, Evaluator > >
    simpleGridFunction ( const GridView &gridView, Evaluator evaluator, int order = std::numeric_limits< int >::max() )
    {
      return SimpleGlobalGridFunction< GridView, Evaluator >( gridView, std::move( evaluator ), order );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_SIMPLEGRIDFUNCTION_HH
