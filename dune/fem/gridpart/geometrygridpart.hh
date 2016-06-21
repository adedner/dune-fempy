#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_HH

#include <dune/common/version.hh>

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>

#if ! DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
#error "Experimental grid extensions required for GeometryGridPart. Add -DDUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS=TRUE to your CMAKE_FLAGS."
#else // #if ! DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS

#include <dune/fem/gridpart/common/deaditerator.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/gridpart/common/metatwistutility.hh>
#include <dune/fem/gridpart/geometrygridpart/capabilities.hh>
#include <dune/fem/gridpart/geometrygridpart/datahandle.hh>
#include <dune/fem/gridpart/geometrygridpart/entity.hh>
#include <dune/fem/gridpart/geometrygridpart/intersection.hh>
#include <dune/fem/gridpart/geometrygridpart/intersectioniterator.hh>
#include <dune/fem/gridpart/geometrygridpart/iterator.hh>
#include <dune/fem/gridpart/geometrygridpart/geometry.hh>
#include <dune/fem/gridpart/idgridpart/indexset.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class GridFunctionType >
    class GeometryGridPart;



    // GeometryGridPartFamily
    // ----------------------

    template< class GridFunction >
    struct GeometryGridPartFamily
    {
      typedef GridFunction GridFunctionType;
      typedef typename GridFunction::RangeFieldType ctype;

      static const int dimension = GridFunction::GridPartType::dimension;
      static const int dimensionworld = GridFunction::FunctionSpaceType::dimRange;

      typedef GeometryGridPartFamily< GridFunction > GridPartFamily;

      struct Traits
      {
        typedef GridFunction GridFunctionType;
        typedef typename GridFunctionType::GridPartType HostGridPartType;

        template< int codim >
        struct Codim
        {
          typedef Dune::Geometry< dimension - codim, dimensionworld, const GridPartFamily, GeometryGridPartGeometry > Geometry;
          typedef typename HostGridPartType::template Codim< codim >::LocalGeometryType LocalGeometry;

          typedef Dune::Entity< codim, dimension, const GridPartFamily, GeometryGridPartEntity > Entity;
          typedef typename HostGridPartType::GridType::template Codim< codim >::EntitySeed EntitySeed;
#if ! DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )
          typedef Dune::EntityPointer< const GridPartFamily, DefaultEntityPointer< Entity > > EntityPointer;
#endif // #if ! DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )
        };

        typedef DeadIntersection< const GridPartFamily > IntersectionImplType;
        typedef DeadIntersectionIterator< const GridPartFamily > IntersectionIteratorImplType;

        typedef Dune::Intersection< const GridPartFamily, IntersectionImplType > LeafIntersection;
        typedef Dune::Intersection< const GridPartFamily, IntersectionImplType > LevelIntersection;

        typedef Dune::IntersectionIterator< const GridPartFamily, IntersectionIteratorImplType, IntersectionImplType > LeafIntersectionIterator;
        typedef Dune::IntersectionIterator< const GridPartFamily, IntersectionIteratorImplType, IntersectionImplType > LevelIntersectionIterator;

        typedef Dune::EntityIterator< 0, const GridPartFamily, DeadIterator< typename Codim< 0 >::Entity > > HierarchicIterator;
      };

      template< int codim >
      struct Codim
      : public Traits::template Codim< codim >
      {};

      typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

      typedef typename Traits::HierarchicIterator HierarchicIterator;
    };

    // GeometryGridPartTraits
    // ----------------------

    template< class GridFunction >
    struct GeometryGridPartTraits
    {
      typedef GridFunction GridFunctionType;
      typedef typename GridFunction :: GridPartType HostGridPart;
      typedef GeometryGridPart< GridFunction > GridPartType;
      typedef GeometryGridPartFamily< GridFunction > GridPartFamily;
      typedef GeometryGridPartFamily< GridFunction > GridFamily;

      typedef GridPart2GridViewImpl< GridPartType > GridViewType;

      static const int dimension = GridFunction::GridPartType::dimension;
      static const int dimensionworld = GridFunction::FunctionSpaceType::dimRange;

      //! type of twist utility
      typedef MetaTwistUtility< typename HostGridPart :: TwistUtilityType >  TwistUtilityType;

      typedef IdIndexSet< const GridPartFamily > IndexSetType;

      typedef typename HostGridPart::GridType GridType;

      static const PartitionIteratorType indexSetPartitionType = HostGridPart::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = HostGridPart::indexSetInterfaceType;

      typedef GeometryGridPartIntersectionIterator < const GridFamily > IntersectionIteratorImplType;
      typedef GeometryGridPartIntersection< const GridFamily > IntersectionImplType;
      typedef IntersectionIterator< const GridFamily, IntersectionIteratorImplType, IntersectionImplType > IntersectionIteratorType;

      template< int codim >
      struct Codim
      {
        typedef typename GridFamily::Traits::template Codim< codim >::Geometry GeometryType;
        typedef typename GridFamily::Traits::template Codim< codim >::LocalGeometry LocalGeometryType;

        typedef typename GridFamily::Traits::template Codim< codim >::Entity EntityType;
        typedef Dune::EntityPointer< const GridPartFamily, DefaultEntityPointer< EntityType > > EntityPointerType;
        typedef typename GridFamily::Traits::template Codim< codim >::EntitySeed EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef EntityIterator< codim, const GridFamily, GeometryGridPartIterator< codim, pitype, const GridFamily > > IteratorType;
        };
      };

      typedef typename HostGridPart::CollectiveCommunicationType CollectiveCommunicationType;
      static const bool conforming = HostGridPart::Traits::conforming;
    };



    // GeometryGridPart
    // ----------

    template< class GridFunction >
    struct GeometryGridPart
    : public GridPartInterface< GeometryGridPartTraits< GridFunction > >
    {
      typedef GridFunction GridFunctionType;
    private:
      typedef GeometryGridPart< GridFunctionType > ThisType;
      typedef GridPartInterface< GeometryGridPartTraits< GridFunctionType > > BaseType;
      typedef typename GeometryGridPartTraits< GridFunctionType >::GridFamily GridFamily;
    public:
      typedef typename GridFunctionType::GridPartType HostGridPart;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::CollectiveCommunicationType CollectiveCommunicationType;

      // the interface takes this from the grid
      static const int dimensionworld = GridFunction::FunctionSpaceType::dimRange;

      template< int codim >
      struct Codim
      : public BaseType::template Codim< codim >
      {};

      explicit GeometryGridPart ( const GridFunctionType &gridFunction )
      : gridFunction_( gridFunction ),
        indexSet_( hostGridPart().indexSet() )
      {}

      const GridType &grid () const
      {
        return hostGridPart().grid();
      }

      GridType &grid ()
      {
        return const_cast< GridType & >( hostGridPart().grid() ); //! correct?
      }

      const IndexSetType &indexSet () const
      {
        return indexSet_;
      }

      template< int codim >
      typename Codim< codim >::IteratorType
      begin () const
      {
        return begin< codim, InteriorBorder_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      begin () const
      {
        return GeometryGridPartIterator< codim, pitype, const GridFamily >( gridFunction_, hostGridPart().template begin< codim, pitype >() );
      }

      template< int codim >
      typename Codim< codim >::IteratorType
      end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      end () const
      {
        return GeometryGridPartIterator< codim, pitype, const GridFamily >( gridFunction_, hostGridPart().template end< codim, pitype >() );
      }

      int level () const
      {
        return hostGridPart().level();
      }

      IntersectionIteratorType ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return GeometryGridPartIntersectionIterator< const GridFamily >( entity, hostGridPart().ibegin( entity.impl().hostEntity() ) );
      }

      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return GeometryGridPartIntersectionIterator< const GridFamily >( entity, hostGridPart().iend( entity.impl().hostEntity() ) );
      }

      int boundaryId ( const IntersectionType &intersection ) const
      {
        return hostGridPart().boundaryId( intersection.impl().hostIntersection() );
      }

      const CollectiveCommunicationType &comm () const { return hostGridPart().comm(); }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &handle,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        typedef CommDataHandleIF< DataHandle, Data >  HostHandleType;
        GeometryGridPartDataHandle< GridFamily, HostHandleType > handleWrapper( handle, gridFunction_ );
        hostGridPart().communicate( handleWrapper, iftype, dir );
      }

#if ! DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )
      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityPointerType
      entityPointer ( const EntitySeed &seed ) const
      {
        typedef typename Codim< EntitySeed::codimension >::EntityPointerType::Implementation EntityPointerImp;
        return EntityPointerImp( hostGridPart().entityPointer( seed ), gridFunction_ );
      }
#endif // #if ! DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )

      // convert a grid entity to a grid part entity ("Gurke!")
      template< class Entity >
      MakeableInterfaceObject< typename Codim< Entity::codimension >::EntityType >
      convert ( const Entity &entity ) const
      {
        // create a grid part entity from a given grid entity
        typedef typename Codim< Entity::codimension >::EntityType EntityType;
        typedef typename EntityType::Implementation Implementation;
        typedef MakeableInterfaceObject< EntityType > EntityObj;
        // here, grid part information can be passed, if necessary
        return EntityObj( Implementation( entity, gridFunction_ ) );
      }

    private:
      const HostGridPart &hostGridPart () const
      {
        return gridFunction_.gridPart();
      }

      const GridFunctionType &gridFunction_;
      IndexSetType indexSet_;
    };



    // GridEntityAccess for GeometryGridPartEntity
    // -------------------------------------------

    template< int codim, int dim, class GridFamily >
    struct GridEntityAccess< Dune::Entity< codim, dim, GridFamily, GeometryGridPartEntity > >
    {
      typedef Dune::Entity< codim, dim, GridFamily, GeometryGridPartEntity > EntityType;
      typedef GridEntityAccess< typename EntityType::Implementation::HostEntityType > HostAccessType;
      typedef typename HostAccessType::GridEntityType GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return HostAccessType::gridEntity( entity.impl().hostEntity() );
      }
    };



    // EntitySearch for GeometryGridPart
    // ---------------------------------

    template< class GridFunction, int codim, PartitionIteratorType partition >
    class EntitySearch< GeometryGridPart< GridFunction >, codim, partition >
      : public DefaultEntitySearch< GeometryGridPart< GridFunction >, codim, partition >
    {
      typedef EntitySearch< GeometryGridPart< GridFunction >, codim, partition > ThisType;
      typedef DefaultEntitySearch< GeometryGridPart< GridFunction >, codim, partition > BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      explicit EntitySearch ( const GridPartType &gridPart )
        : BaseType( gridPart )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #else // #if ! DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_HH