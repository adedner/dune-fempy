#ifndef DUNE_FEMPY_GRID_HIERARCHICAL_HH
#define DUNE_FEMPY_GRID_HIERARCHICAL_HH

#include <cstddef>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>

#include <dune/fempy/grid/restrictprolong.hh>
#include <dune/fempy/parameter.hh>

namespace Dune
{

  namespace FemPy
  {

    /*!
        @file
        @brief Contains C++ template classes for grids.
        \ingroup Grids
    */

    // readDGF
    // -------

    template< class Grid >
    static std::shared_ptr< Grid > readDGF ( const std::string &dgf )
    {
      GridPtr< Grid > gridPtr( dgf );
      gridPtr->loadBalance();
      return std::shared_ptr< Grid >( gridPtr.release() );
    }



    // DiscreteFunctionList
    // --------------------

    template< class Grid, class D = double >
    struct DiscreteFunctionList
      : public Fem::IsDiscreteFunction
    {
      typedef AdaptiveDofVector< Grid, D > DiscreteFunction;

      typedef typename DiscreteFunction::DofType DofType;

      typedef typename Grid::template Codim< 0 >::Entity ElementType;

      struct GridPartType
      {
        explicit GridPartType ( std::shared_ptr< Grid > grid ) : grid_( std::move( grid ) ) {}

        struct IndexSetType
        {
          bool contains ( const ElementType & ) const { return true; }
        };

        typedef Grid GridType;

        const Grid &grid () const { return *grid_; }

        const ElementType &convert ( const ElementType &entity ) const { return entity; }

      private:
        std::shared_ptr< Grid > grid_;
      };

      struct DiscreteFunctionSpaceType
      {
        typedef std::vector< std::shared_ptr< DiscreteFunction > > DiscreteFunctions;

        typedef DiscreteFunctionList::GridPartType GridPartType;

        typedef typename GridPartType::IndexSetType IndexSetType;
        typedef typename GridPartType::GridType GridType;
        typedef ElementType EntityType;

        struct BasisFunctionSetType
        {
          BasisFunctionSetType ( const DiscreteFunctions &discreteFunctions, const ElementType &entity )
            : discreteFunctions_( discreteFunctions ), entity_( entity )
          {}

          std::size_t size () const
          {
            std::size_t size = 0;
            for( const auto &df : discreteFunctions_ )
              size += df->numLocalDofs( entity_ );
            return size;
          }

        private:
          const DiscreteFunctions &discreteFunctions_;
          const ElementType &entity_;
        };

        explicit DiscreteFunctionSpaceType ( GridPartType gridPart ) : gridPart_( std::move( gridPart ) ) {}

        const GridPartType &gridPart () const { return gridPart_; }
        IndexSetType indexSet () const { return IndexSetType(); }

        BasisFunctionSetType basisFunctionSet ( const ElementType &entity ) const { return BasisFunctionSetType( discreteFunctions_, entity ); }

        const DiscreteFunctions &discreteFunctions () const { return discreteFunctions_; }
        DiscreteFunctions &discreteFunctions () { return discreteFunctions_; }

      private:
        GridPartType gridPart_;
        DiscreteFunctions discreteFunctions_;
      };

      typedef typename DiscreteFunctionSpaceType::DiscreteFunctions::const_iterator ConstIterator;
      typedef typename DiscreteFunctionSpaceType::DiscreteFunctions::iterator Iterator;

      typedef std::allocator< DofType > LocalDofVectorAllocatorType;

      explicit DiscreteFunctionList ( std::shared_ptr< Grid > grid )
        : space_( GridPartType( std::move( grid ) ) )
      {}

      const GridPartType &gridPart () const { return space_.gridPart(); }

      const DiscreteFunctionSpaceType &space () const { return space_; }

      void enableDofCompression ()
      {
        for( const auto &df : space_.discreteFunctions() )
          df->enableDofCompression();
      }

      void assign () { space_.discreteFunctions().clear(); }

      template< class Iterator >
      void assign ( Iterator begin, Iterator end )
      {
        space_.discreteFunctions().assign( begin, end );
      }

      ConstIterator begin () const { return space_.discreteFunctions().begin(); }
      Iterator begin () { return space_.discreteFunctions().begin(); }
      ConstIterator end () const { return space_.discreteFunctions().end(); }
      Iterator end () { return space_.discreteFunctions().end(); }

      LocalDofVectorAllocatorType localDofVectorAllocator () const { return LocalDofVectorAllocatorType(); }

      template< class A >
      void getLocalDofs ( const ElementType &entity, DynamicVector< DofType, A > &localDofs ) const
      {
        DofType *it = &localDofs[ 0 ];
        for( const auto &df : space_.discreteFunctions() )
          it = df->getLocalDofs( entity, it );
      }

      template< class A >
      void setLocalDofs ( const ElementType &entity, const DynamicVector< DofType, A > &localDofs )
      {
        const DofType *it = &localDofs[ 0 ];
        for( auto &df : space_.discreteFunctions() )
          it = df->setLocalDofs( entity, it );
      }

    private:
      DiscreteFunctionSpaceType space_;
    };

  } // namespace FemPy



  namespace Fem
  {

    // DiscreteFunctionTraits for DiscreteFunctionList
    // -----------------------------------------------

    template< class Grid, class D >
    struct DiscreteFunctionTraits< FemPy::DiscreteFunctionList< Grid, D > >
    {
      typedef typename FemPy::AdaptiveDofVector< Grid, D >::DofType DofType;
      typedef std::allocator< DofType > LocalDofVectorAllocatorType;
    };

  } // namespace Fem



  namespace FemPy
  {

    // RestrictProlong
    // ---------------

    template< class Grid, class D = double >
    class RestrictProlong
      : public Fem::RestrictProlongInterface< Fem::RestrictProlongTraits< RestrictProlong< Grid, D >, typename Grid::ctype > >
    {
      typedef Fem::RestrictProlongInterface< Fem::RestrictProlongTraits< RestrictProlong< Grid, D >, typename Grid::ctype > > BaseType;

      struct LoadBalanceContainsCheck;

    public:
      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef typename Grid::template Codim< 0 >::Entity ElementType;
      typedef typename Grid::template Codim< 0 >::LocalGeometry LocalGeometryType;

      explicit RestrictProlong ( std::shared_ptr< Grid > grid )
        : discreteFunctions_( std::move( grid ) )
      {}

      void setFatherChildWeight ( const DomainFieldType &weight ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().setFatherChildWeight( weight );
      }

      void restrictLocal ( const ElementType &father, const ElementType &child, bool initialize ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().restrictLocal( father, child, initialize );
      }

      void restrictLocal ( const ElementType &father, const ElementType &child, const LocalGeometryType &geometryInFather, bool initialize ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().restrictLocal( father, child, geometryInFather, initialize );
      }

      void prolongLocal ( const ElementType &father, const ElementType &child, bool initialize ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().prolongLocal( father, child, initialize );
      }

      void prolongLocal ( const ElementType &father, const ElementType &child, const LocalGeometryType &geometryInFather, bool initialize ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().prolongLocal( father, child, geometryInFather, initialize );
      }

      void addToList ( Fem::CommunicationManagerList &commList )
      {
        if( commList_ )
          DUNE_THROW( InvalidStateException, "Only one communication list supported." );
        commList_ = &commList;
        addToCommList();
      }

      void removeFromList ( Fem::CommunicationManagerList &commList )
      {
        if( commList_ != &commList )
          DUNE_THROW( InvalidStateException, "Only one communication list supported." );
        removeFromCommList();
        commList_ = nullptr;
      }

      template< class LoadBalancer >
      void addToLoadBalancer ( LoadBalancer &loadBalancer );

      void assign ()
      {
        removeFromCommList();
        discreteFunctions_.assign();
        // no need to add to comm list; discrete functions is empty
      }

      template< class Iterator >
      void assign ( Iterator begin, Iterator end )
      {
        removeFromCommList();
        discreteFunctions_.assign( std::move( begin ), std::move( end ) );
        addToCommList();
      }

    private:
      void addToCommList ()
      {
        if( commList_ )
          for( const auto &df : discreteFunctions_ )
            df->restrictProlong().addToList( *commList_ );
      }

      void removeFromCommList ()
      {
        if( commList_ )
          for( const auto &df : discreteFunctions_ )
            df->restrictProlong().removeFromList( *commList_ );
      }

      DiscreteFunctionList< Grid, D > discreteFunctions_;
      Fem::CommunicationManagerList *commList_ = nullptr;
    };



    // RestrictProlong::LoadBalanceContainsCheck
    // -----------------------------------------

    template< class Grid, class D >
    struct RestrictProlong< Grid, D >::LoadBalanceContainsCheck
    {
      explicit LoadBalanceContainsCheck ( const DiscreteFunctionList< Grid, D > &discreteFunctions )
        : discreteFunctions_( discreteFunctions )
      {}

      bool contains ( const ElementType &element ) const
      {
        // todo: implement this predicate
        return false;
      }

    private:
      const DiscreteFunctionList< Grid, D > &discreteFunctions_;
    };



    // Implementation of RestrictProlong
    // ---------------------------------

    template< class Grid, class D >
    template< class LoadBalancer >
    inline void RestrictProlong< Grid, D >::addToLoadBalancer ( LoadBalancer &loadBalancer )
    {
      loadBalancer.addDiscreteFunction( discreteFunctions_, LoadBalanceContainsCheck( discreteFunctions_ ) );
    }



    // HierarchicalGrid
    // ----------------

    template< class G >
    struct HierarchicalGrid
    {
      typedef G Grid;

      typedef Fem::AdaptationManager< Grid, RestrictProlong< Grid > > AdaptationManager;

      enum class Marker { Coarsen = -1, Keep = 0, Refine = 1 };

      typedef typename Grid::template Codim< 0 >::Entity Element;

      explicit HierarchicalGrid ( const std::string &dgf )
        : grid_( readDGF< Grid >( dgf ) ),
          restrictProlong_( new RestrictProlong< Grid >( grid_ ) ),
          adaptationManager_( new AdaptationManager( *grid_, *restrictProlong_, noParameter() ) )
      {}

      template< class Marking >
      void mark ( Marking marking )
      {
        for( const Element &element : elements( grid_->leafGridView() ) )
        {
          Marker marker = marking( element );
          grid_->mark( static_cast< int >( marker ), element );
        }
      }

      template< class Iterator >
      void adapt ( Iterator begin, Iterator end )
      {
        restrictProlong_->assign( begin, end );
        adaptationManager_->adapt();
        restrictProlong_->assign();
      }

      void globalRefine ( int level )
      {
        Fem::GlobalRefine::apply( *grid_, level );
      }

      template< class Iterator >
      void loadBalance ( Iterator begin, Iterator end )
      {
        restrictProlong_->assign( begin, end );
        adaptationManager_->loadBalance();
        restrictProlong_->assign();
      }

      std::shared_ptr< Grid > grid () const { return grid_; }

    private:
      std::shared_ptr< Grid > grid_;
      std::shared_ptr< RestrictProlong< Grid > > restrictProlong_;
      std::shared_ptr< AdaptationManager > adaptationManager_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_HIERARCHICAL_HH