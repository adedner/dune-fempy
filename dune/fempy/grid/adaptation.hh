#ifndef DUNE_FEMPY_GRID_ADAPTATION_HH
#define DUNE_FEMPY_GRID_ADAPTATION_HH

#include <cstddef>

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>

#include <dune/fempy/grid/adaptivedofvector.hh>
#include <dune/fempy/grid/virtualizedrestrictprolong.hh>
#include <dune/fempy/parameter.hh>

namespace Dune
{

  namespace FemPy
  {

    // FakeGridPart
    // ------------

    template< class Grid >
    class FakeGridPart
    {
      typedef typename Grid::template Codim< 0 >::Entity Element;

    public:
      explicit FakeGridPart ( Grid &grid ) : grid_( grid ) {}

      struct IndexSetType
      {
        bool contains ( const Element & ) const { return true; }
      };

      typedef Grid GridType;

      Grid &grid () const { return grid_; }

      const Element &convert ( const Element &entity ) const { return entity; }

    private:
      Grid &grid_;
    };



    // FakeDiscreteFunctionSpace
    // -------------------------

    template< class Grid, class NumLocalDofs = std::function< std::size_t( typename Grid::template Codim< 0 >::Entity ) > >
    struct FakeDiscreteFunctionSpace
    {
      typedef FakeGridPart< Grid > GridPartType;

      typedef typename GridPartType::IndexSetType IndexSetType;
      typedef typename GridPartType::GridType GridType;

      typedef typename GridType::template Codim< 0 >::Entity EntityType;

      struct BasisFunctionSetType
      {
        BasisFunctionSetType ( NumLocalDofs numLocalDofs, const EntityType &entity ) : numLocalDofs_( std::move( numLocalDofs ) ), entity_( entity ) {}

        std::size_t size () const { return numLocalDofs_( entity_ ); }

      private:
        NumLocalDofs numLocalDofs_;
        const EntityType &entity_;
      };

      explicit FakeDiscreteFunctionSpace ( GridType &grid, NumLocalDofs numLocalDofs ) : gridPart_( grid ), numLocalDofs_( std::move( numLocalDofs ) ) {}
      explicit FakeDiscreteFunctionSpace ( GridPartType gridPart, NumLocalDofs numLocalDofs ) : gridPart_( std::move( gridPart ) ), numLocalDofs_( std::move( numLocalDofs ) ) {}

      const GridPartType &gridPart () const { return gridPart_; }
      IndexSetType indexSet () const { return IndexSetType(); }

      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const { return BasisFunctionSetType( numLocalDofs_, entity ); }

    private:
      GridPartType gridPart_;
      NumLocalDofs numLocalDofs_;
    };



    // fakeDiscreteFunctionSpace
    // -------------------------

    template< class Grid, class NumLocalDofs >
    inline static FakeDiscreteFunctionSpace< Grid, NumLocalDofs > fakeDiscreteFunctionSpace ( const Grid &grid, NumLocalDofs numLocalDofs )
    {
      return FakeDiscreteFunctionSpace< Grid, NumLocalDofs >( grid, std::move( numLocalDofs ) );
    }



    // DiscreteFunctionList
    // --------------------

    template< class Grid, class D = double >
    struct DiscreteFunctionList
      : public Fem::IsDiscreteFunction
    {
      typedef AdaptiveDofVector< Grid, D > DiscreteFunction;

      typedef typename DiscreteFunction::DofType DofType;

      typedef FakeDiscreteFunctionSpace< Grid > DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::EntityType ElementType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      typedef std::vector< std::shared_ptr< DiscreteFunction > > DiscreteFunctions;

      typedef typename DiscreteFunctions::const_iterator ConstIterator;
      typedef typename DiscreteFunctions::iterator Iterator;

      typedef std::allocator< DofType > LocalDofVectorAllocatorType;

      explicit DiscreteFunctionList ( Grid &grid )
        : space_( grid, [ this ] ( const ElementType &element ) {
              std::size_t size = 0;
              for( const auto &df : discreteFunctions_ )
                size += df->numLocalDofs( element );
              return size;
            } )
      {}

      const GridPartType &gridPart () const { return space_.gridPart(); }

      const DiscreteFunctionSpaceType &space () const { return space_; }

      void enableDofCompression ()
      {
        for( const auto &df : discreteFunctions_ )
          df->enableDofCompression();
      }

      void assign () { discreteFunctions_.clear(); }

      template< class Iterator >
      void assign ( Iterator begin, Iterator end )
      {
        discreteFunctions_.assign( begin, end );
      }

      ConstIterator begin () const { return discreteFunctions_.begin(); }
      Iterator begin () { return discreteFunctions_.begin(); }
      ConstIterator end () const { return discreteFunctions_.end(); }
      Iterator end () { return discreteFunctions_.end(); }

      LocalDofVectorAllocatorType localDofVectorAllocator () const { return LocalDofVectorAllocatorType(); }

      template< class A >
      void getLocalDofs ( const ElementType &entity, DynamicVector< DofType, A > &localDofs ) const
      {
        DofType *it = &localDofs[ 0 ];
        for( const auto &df : discreteFunctions_ )
          it = df->getLocalDofs( entity, it );
      }

      template< class A >
      void setLocalDofs ( const ElementType &entity, const DynamicVector< DofType, A > &localDofs )
      {
        const DofType *it = &localDofs[ 0 ];
        for( auto &df : discreteFunctions_ )
          it = df->setLocalDofs( entity, it );
      }

    private:
      DiscreteFunctionSpaceType space_;
      DiscreteFunctions discreteFunctions_;
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

      explicit RestrictProlong ( Grid &grid )
        : discreteFunctions_( grid )
      {}

      void setFatherChildWeight ( const DomainFieldType &weight ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.setFatherChildWeight( weight );
      }

      void restrictLocal ( const ElementType &father, const ElementType &child, bool initialize ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.restrictLocal( father, child, initialize );
      }

      void restrictLocal ( const ElementType &father, const ElementType &child, const LocalGeometryType &geometryInFather, bool initialize ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.restrictLocal( father, child, geometryInFather, initialize );
      }

      void prolongLocal ( const ElementType &father, const ElementType &child, bool initialize ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.prolongLocal( father, child, initialize );
      }

      void prolongLocal ( const ElementType &father, const ElementType &child, const LocalGeometryType &geometryInFather, bool initialize ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.prolongLocal( father, child, geometryInFather, initialize );
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
        restrictProlongs_.clear();
        discreteFunctions_.assign();
        // no need to add to comm list; discrete functions is empty
      }

      template< class Iterator >
      void assign ( Iterator begin, Iterator end )
      {
        removeFromCommList();
        restrictProlongs_.clear();
        discreteFunctions_.assign( std::move( begin ), std::move( end ) );
        for( const auto &df : discreteFunctions_ )
          restrictProlongs_.push_back( df->restrictProlong() );
        addToCommList();
      }

      Grid &grid () const { return discreteFunctions_.gridPart().grid(); }

    private:
      void addToCommList ()
      {
        if( commList_ )
          for( auto &rp : restrictProlongs_ )
            rp.addToList( *commList_ );
      }

      void removeFromCommList ()
      {
        if( commList_ )
          for( auto &rp : restrictProlongs_ )
            rp.removeFromList( *commList_ );
      }

      DiscreteFunctionList< Grid, D > discreteFunctions_;
      std::vector< VirtualizedRestrictProlong< Grid > > restrictProlongs_;
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



    // GridAdaptation
    // --------------

    template< class G >
    struct GridAdaptation
    {
      typedef G Grid;

      typedef Fem::AdaptationManager< Grid, RestrictProlong< Grid > > AdaptationManager;

      enum class Marker { Coarsen = -1, Keep = 0, Refine = 1 };

      typedef typename Grid::template Codim< 0 >::Entity Element;

      explicit GridAdaptation ( Grid &grid )
        : restrictProlong_( grid ),
          adaptationManager_( grid, restrictProlong_, noParameter() )
      {}

      template< class Marking >
      void mark ( Marking marking )
      {
        for( const Element &element : elements( grid().leafGridView() ) )
        {
          Marker marker = marking( element );
          grid().mark( static_cast< int >( marker ), element );
        }
      }

      template< class Iterator >
      void adapt ( Iterator begin, Iterator end )
      {
        restrictProlong_.assign( begin, end );
        adaptationManager_.adapt();
        restrictProlong_.assign();
      }

      void globalRefine ( int level )
      {
        Fem::GlobalRefine::apply( grid(), level );
      }

      template< class Iterator >
      void loadBalance ( Iterator begin, Iterator end )
      {
        restrictProlong_.assign( begin, end );
        adaptationManager_.loadBalance();
        restrictProlong_.assign();
      }

      Grid &grid () const { return restrictProlong_.grid(); }

    private:
      RestrictProlong< Grid > restrictProlong_;
      AdaptationManager adaptationManager_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_ADAPTATION_HH