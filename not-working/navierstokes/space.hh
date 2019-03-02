#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_MINIELEMENT_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_MINIELEMENT_HH

#include <dune/grid/common/gridenums.hh>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>

#include <dune/fem/space/basisfunctionset/default.hh>

#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

#include <dune/fem/space/mapper/indexsetdofmapper.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>

#include "localinterpolation.hh"
#include "localkeymap.hh"
#include "shapefunctions.hh"

namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    // -------------------

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    class BubbleElementSpace;

    // BubbleElementSpaceTraits
    // ----------------------

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    struct BubbleElementSpaceTraits
    {
      typedef BubbleElementSpace< FunctionSpace, GridPart, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
    public:
      typedef typename FunctionSpaceType :: ScalarFunctionSpaceType ScalarFunctionSpaceType;

      // defined in shapefunctionset.hh
      typedef SimplexBubbleElementShapeFunctionSet< ScalarFunctionSpaceType > ScalarShapeFunctionSetType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;

      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      static const int localBlockSize = FunctionSpaceType::dimRange;
      static const int polynomialOrder = ScalarShapeFunctionSetType::polynomialOrder;

      typedef IndexSetDofMapper< GridPartType > BlockMapperType;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };

    // BubbleElementSpace
    // ----------------

    //! [Class definition for new space]
    template< class FunctionSpace, class GridPart, template< class > class Storage = CachingStorage >
    class BubbleElementSpace
    : public DiscreteFunctionSpaceDefault< BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > >
    //! [Class definition for new space]
    {
      typedef BubbleElementSpace< FunctionSpace, GridPart, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > > BaseType;

    public:
      typedef BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > Traits;
      typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;
      static const int polynomialOrder = Traits::polynomialOrder;

      typedef typename BaseType::GridPartType GridPartType;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      // type of local interpolation
      typedef LocalBubbleElementInterpolation< FunctionSpace > InterpolationType;

      // static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;
      static const InterfaceType defaultInterface = GridPartType::indexSetInterfaceType;
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      BubbleElementSpace ( GridPartType &gridPart,
                         const InterfaceType commInterface = defaultInterface,
                         const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        blockMapper_(  gridPart, BubbleElementDofMapperCodeFactory() )
      {
      }

      ~BubbleElementSpace ()
      {
      }

      bool contains ( const int codim ) const
      {
        // forward to mapper since this information is held there
        return blockMapper().contains( codim );
      }

      bool continuous () const { return true; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType & intersection ) const
      {
        // forward to the subsapces
        return true;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      template<class Entity>
      int order ( const Entity &entity ) const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet(const EntityType &entity) const */
      template< class EntityType >
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, ShapeFunctionSetType( entity.geometry().type() ) );
      }

      /** \brief obtain the DoF block mapper of this space
          \return BlockMapperType
      **/
      BlockMapperType &blockMapper () const { return blockMapper_; }

      // Non-interface methods

      /** \brief return local interpolation for given entity
       *
       *  \param[in]  entity  grid part entity
       *  \note this method is needed to call the global function inteprolate( ... )
       */
      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return InterpolationType();
      }

    protected:
      mutable BlockMapperType blockMapper_;
    };

    template< class FunctionSpace,
              class GridPart,
              template<class> class Storage,
              class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< BubbleElementSpace < FunctionSpace, GridPart, Storage >, NewFunctionSpace >
    {
      typedef BubbleElementSpace< NewFunctionSpace, GridPart, Storage > Type;
    };

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    class DefaultLocalRestrictProlong < BubbleElementSpace< FunctionSpace, GridPart, Storage > >
    : public EmptyLocalRestrictProlong< BubbleElementSpace< FunctionSpace, GridPart, Storage > >
    {
      typedef EmptyLocalRestrictProlong< BubbleElementSpace< FunctionSpace, GridPart, Storage > > BaseType;
      public:
      DefaultLocalRestrictProlong( const BubbleElementSpace< FunctionSpace, GridPart, Storage > &space )
        : BaseType()
      {}
    };

    template< class MatrixImp,
              class FunctionSpace, class GridPart, template< class > class Storage>
    struct ISTLParallelMatrixAdapter< MatrixImp, BubbleElementSpace< FunctionSpace,GridPart,Storage> >
    {
      typedef LagrangeParallelMatrixAdapter<MatrixImp> Type;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_MINIELEMENT_HH
