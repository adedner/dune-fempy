/**************************************************************************

  The dune-fem module is a module of DUNE (see www.dune-project.org).
  It is based on the dune-grid interface library
  extending the grid interface by a number of discretization algorithms
  for solving non-linear systems of partial differential equations.

  Copyright (C) 2003 - 2015 Robert Kloefkorn
  Copyright (C) 2003 - 2010 Mario Ohlberger
  Copyright (C) 2004 - 2015 Andreas Dedner
  Copyright (C) 2005        Adrian Burri
  Copyright (C) 2005 - 2015 Mirko Kraenkel
  Copyright (C) 2006 - 2015 Christoph Gersbacher
  Copyright (C) 2006 - 2015 Martin Nolte
  Copyright (C) 2011 - 2015 Tobias Malkmus
  Copyright (C) 2012 - 2015 Stefan Girke
  Copyright (C) 2013 - 2015 Claus-Justus Heine
  Copyright (C) 2013 - 2014 Janick Gerstenberger
  Copyright (C) 2013        Sven Kaulman
  Copyright (C) 2013        Tom Ranner
  Copyright (C) 2015        Marco Agnese
  Copyright (C) 2015        Martin Alkaemper


  The dune-fem module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The dune-fem module is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

**************************************************************************/
#ifndef DUNE_DIRICHLETCONSTRAINTS_HH
#define DUNE_DIRICHLETCONSTRAINTS_HH

#include <dune/fem/function/common/scalarproducts.hh>

namespace Dune {

template < class Model, class DiscreteFunctionSpace >
class DirichletConstraints
{
public:
  typedef Model ModelType;
  typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

  //! type of grid partition
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  //! type of grid
  typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  // types for boundary treatment
  // ----------------------------
  typedef typename DiscreteFunctionSpaceType :: BlockMapperType BlockMapperType;

  typedef Fem::SlaveDofs< DiscreteFunctionSpaceType, BlockMapperType > SlaveDofsType;
  typedef typename SlaveDofsType :: SingletonKey SlaveDofsKeyType;
  typedef Fem::SingletonList< SlaveDofsKeyType, SlaveDofsType >
      SlaveDofsProviderType;

  static const int localBlockSize = DiscreteFunctionSpaceType :: localBlockSize ;
  static_assert( localBlockSize == DiscreteFunctionSpaceType::FunctionSpaceType::dimRange,
      "local block size of the space must be identical to the dimension of the range of the function space.");
  typedef FieldVector<int, localBlockSize> DirichletBlock;
  typedef FieldVector<int, ModelType::dimRange> ModelDirichletBlock;
  static_assert( ModelType::dimRange >= localBlockSize,
      "local block size of the space must be less or equahl to the dimension of the range of the model.");

  DirichletConstraints( const ModelType &model, const DiscreteFunctionSpaceType& space )
    : model_(model),
      space_( space ),
      slaveDofs_( getSlaveDofs( space_ ) ),
      dirichletBlocks_(),
      // mark DoFs on the Dirichlet boundary
      hasDirichletDofs_( false ),
      sequence_( -1 )
  {
  }

  /*! treatment of Dirichlet-DoFs for given discrete function
   *
   *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints
   *   \param[out] w   discrete function the constraints are applied to
   */
  template < class DiscreteFunctionType >
  void operator ()( const DiscreteFunctionType& u, DiscreteFunctionType& w ) const
  {
    updateDirichletDofs();

    // if Dirichlet Dofs have been found, treat them
    if( hasDirichletDofs_ )
    {
      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType ;
      typedef typename DiscreteFunctionType :: ConstDofIteratorType ConstDofIteratorType ;

      ConstDofIteratorType uIt = u.dbegin();
      DofIteratorType wIt = w.dbegin();

      // loop over all blocks
      const unsigned int blocks = space_.blockMapper().size();
      for( unsigned int blockDof = 0; blockDof < blocks ; ++ blockDof )
      {
        for( int l = 0; l < localBlockSize ; ++ l, ++ wIt, ++ uIt )
        {
          if( dirichletBlocks_[ blockDof ][l] )
          {
            // copy dofs of the block
            assert( uIt != u.dend() );
            assert( wIt != w.dend() );
            (*wIt) = (*uIt);
          }
        }
      }
    }
  }

  /*! treatment of Dirichlet-DoFs for given discrete function
   *
   *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints
   *   \param[out] w   discrete function the constraints are applied to
   */
  template < class DiscreteFunctionType >
  void operator ()( DiscreteFunctionType& w ) const
  {
    updateDirichletDofs();

    if( hasDirichletDofs_ )
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename IteratorType :: Entity EntityType;

      for( const EntityType &entity : space_ )
      {
        typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

        model_.init(entity);
        LocalFunctionType wLocal = w.localFunction( entity );

        // interpolate dirichlet dofs
        dirichletDofTreatment( wLocal );
      }
    }
  }
  template < class GridFunctionType, class DiscreteFunctionType >
  void operator ()( const GridFunctionType &u, DiscreteFunctionType& w ) const
  {
    updateDirichletDofs();

    if( hasDirichletDofs_ )
    {
      for( const auto &entity : space_ )
      {
        typedef typename GridFunctionType :: LocalFunctionType GridLocalFunctionType;
        typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
        LocalFunctionType wLocal = w.localFunction( entity );
        const GridLocalFunctionType uLocal = u.localFunction( entity );
        // interpolate dirichlet dofs
        dirichletDofTreatment( uLocal, wLocal );
      }
    }
  }

  /*! treatment of Dirichlet-DoFs for solution and right-hand-side
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
   *
   *   \param[out] linearOperator  linear operator to be adjusted
   */
  template <class LinearOperator>
  void applyToOperator( LinearOperator& linearOperator ) const
  {
    updateDirichletDofs();

    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;

    // if Dirichlet Dofs have been found, treat them
    if( hasDirichletDofs_ )
    {
      const IteratorType end = space_.end();
      for( IteratorType it = space_.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;
        // adjust linear operator
        dirichletDofsCorrectOnEntity( linearOperator, entity );
      }
    }
  }

protected:

  /*! treatment of Dirichlet-DoFs for one entity
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
   *
   *   \param[in]  entity  entity to perform Dirichlet treatment on
   */
  template< class LinearOperator, class EntityType >
  void dirichletDofsCorrectOnEntity ( LinearOperator& linearOperator,
                                      const EntityType &entity ) const
  {
    // get slave dof structure (for parallel runs)   /*@LST0S@*/
    SlaveDofsType &slaveDofs = this->slaveDofs();

    typedef typename LinearOperator :: LocalMatrixType LocalMatrixType;

    // get local matrix from linear operator
    LocalMatrixType localMatrix = linearOperator.localMatrix( entity, entity );

    // get number of basis functions
    const int localBlocks = space_.blockMapper().numDofs( entity );

    // map local to global dofs
    std::vector<std::size_t> globalBlockDofs(localBlocks);
    // obtain all DofBlocks for this element
    space_.blockMapper().map( entity, globalBlockDofs );

    // counter for all local dofs (i.e. localBlockDof * localBlockSize + ... )
    int localDof = 0;
    // iterate over face dofs and set unit row
    for( int localBlockDof = 0 ; localBlockDof < localBlocks; ++ localBlockDof )
    {
      int global = globalBlockDofs[localBlockDof];
      for( int l = 0; l < localBlockSize; ++ l, ++ localDof )
      {
        if( dirichletBlocks_[global][l] )
        {
          // clear all other columns
          localMatrix.clearRow( localDof );

          // set diagonal to 1
          double value = slaveDofs.isSlave( global )? 0.0 : 1.0;
          localMatrix.set( localDof, localDof, value );
        }
      }
    }
  }

  //! set the dirichlet points to exact values
  template< class LocalFunctionType >
  void dirichletDofTreatment( LocalFunctionType &wLocal ) const
  {
    // get entity
    const typename LocalFunctionType::EntityType &entity = wLocal.entity();

    // get number of Lagrange Points
    const int localBlocks = space_.blockMapper().numDofs( entity );

    // map local to global BlockDofs
    std::vector<std::size_t> globalBlockDofs(localBlocks);
    space_.blockMapper().map(entity,globalBlockDofs);
    std::vector<typename LocalFunctionType::RangeFieldType> values( localBlocks*localBlockSize );

    int localDof = 0;

    // iterate over face dofs and set unit row
    for( int localBlock = 0 ; localBlock < localBlocks; ++ localBlock )
    {
      // store result to dof vector
      int global = globalBlockDofs[localBlock];
      for( int l = 0; l < localBlockSize ; ++ l, ++localDof )
      {
        if( dirichletBlocks_[ global ][l] )
        {
          space_.interpolation(entity)
            (typename ModelType::BoundaryWrapper(model_, dirichletBlocks_[global][l]), values);
          // store result
          assert( (unsigned int)localDof < wLocal.size() );
          wLocal[ localDof ] = values[ localDof ];
        }
      }
    }
  }
  template< class GridLocalFunctionType, class LocalFunctionType >
  void dirichletDofTreatment( const GridLocalFunctionType &uLocal,
                              LocalFunctionType &wLocal ) const
  {
    // get entity
    const typename LocalFunctionType::EntityType &entity = wLocal.entity();

    // get number of Lagrange Points
    const int localBlocks = space_.blockMapper().numDofs( entity );

    // map local to global BlockDofs
    std::vector<std::size_t> globalBlockDofs(localBlocks);
    space_.blockMapper().map(entity,globalBlockDofs);
    std::vector<double> values( localBlocks*localBlockSize );
    space_.interpolation(entity)(uLocal, values);

    int localDof = 0;

    // iterate over face dofs and set unit row
    for( int localBlock = 0 ; localBlock < localBlocks; ++ localBlock )
    {
      // store result to dof vector
      int global = globalBlockDofs[localBlock];
      for( int l = 0; l < localBlockSize ; ++ l, ++localDof )
      {
        if( dirichletBlocks_[ global ][l] )
        {
          // store result
          assert( (unsigned int)localDof < wLocal.size() );
          wLocal[ localDof ] = values[ localDof ];
        }
      }
    }
  }

protected:
  // detect all DoFs on the Dirichlet boundary
  void updateDirichletDofs() const
  {
    if( sequence_ != space_.sequence() )
    {
      // only start search if Dirichlet boundary is present
      if( ! model_.hasDirichletBoundary() )
      {
        hasDirichletDofs_ = false ;
        return ;
      }

      // resize flag vector with number of blocks and reset flags
      const int blocks = space_.blockMapper().size() ;
      dirichletBlocks_.resize( blocks );
      for( int i=0; i<blocks; ++i )
        dirichletBlocks_[ i ] = DirichletBlock(0) ;

      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename IteratorType :: Entity EntityType;

      bool hasDirichletBoundary = false;
      const IteratorType end = space_.end();
      for( IteratorType it = space_.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;
        // if entity has boundary intersections
        if( entity.hasBoundaryIntersections() )
        {
          hasDirichletBoundary |= searchEntityDirichletDofs( entity, model_ );
        }
      }

      // update sequence number
      sequence_ = space_.sequence();
      if( space_.gridPart().comm().size() > 1 )
      {
        try
        {
          DirichletBuilder handle( *this, space_ , space_.blockMapper() );
          space_.gridPart().communicate
            ( handle, GridPartType::indexSetInterfaceType, ForwardCommunication );
        }
        // catch possible exceptions here to have a clue where it happend
        catch( const Exception &e )
        {
          std::cerr << e << std::endl;
          std::cerr << "Exception thrown in: " << __FILE__ << " line:" << __LINE__ << std::endl;
          abort();
        }
        hasDirichletDofs_ = space_.gridPart().grid().comm().max( hasDirichletBoundary );
      }
      else
      {
        hasDirichletDofs_ = hasDirichletBoundary;
      }
    }
  }

  // detect all DoFs on the Dirichlet boundary of the given entity
  template< class EntityType >
  bool searchEntityDirichletDofs( const EntityType &entity, const ModelType& model ) const
  {
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

    typedef typename GridPartType :: IntersectionIteratorType
      IntersectionIteratorType;

    const GridPartType &gridPart = space_.gridPart();

    // default is false
    bool hasDirichletBoundary = false;

    //map local to global BlockDofs
    std::vector< size_t> globalBlockDofs(space_.blockMapper().numDofs(entity));
    space_.blockMapper().map(entity,globalBlockDofs);

    std::vector< bool>   globalBlockDofsFilter(space_.blockMapper().numDofs(entity));

    IntersectionIteratorType it = gridPart.ibegin( entity );
    const IntersectionIteratorType endit = gridPart.iend( entity );
    for( ; it != endit; ++it )
    {
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;
      const IntersectionType& intersection = *it;

      // if intersection is with boundary, adjust data
      if( intersection.boundary() )
      {
        model.initIntersection( intersection );
        // get dirichlet information from model
        ModelDirichletBlock block(0);
        const bool isDirichletIntersection = model.isDirichletIntersection( intersection, block );
        if (isDirichletIntersection)
        {
          // get face number of boundary intersection
          const int face = intersection.indexInInside();
          space_.blockMapper().onSubEntity(entity,face,1,globalBlockDofsFilter);
          for( unsigned int i=0;i<globalBlockDofs.size();++i)
          {
            if ( !globalBlockDofsFilter[i] ) continue;
            // mark global DoF number
            for(int k = 0; k < DirichletBlock::dimension; ++k)
              dirichletBlocks_[globalBlockDofs[ i ] ][k] = block [k];

            // we have Dirichlet values
            hasDirichletBoundary = true ;
          }
        }
      }
    }

    return hasDirichletBoundary;
  }

  //! pointer to slave dofs
  const ModelType& model_;
  const DiscreteFunctionSpaceType& space_;
  SlaveDofsType *const slaveDofs_;
  mutable std::vector< DirichletBlock > dirichletBlocks_;
  mutable bool hasDirichletDofs_ ;
  mutable int sequence_ ;

  // return slave dofs
  static SlaveDofsType *getSlaveDofs ( const DiscreteFunctionSpaceType &space )
  {
    SlaveDofsKeyType key( space, space.blockMapper() );
    return &(SlaveDofsProviderType :: getObject( key ));
  }

  // return reference to slave dofs
  SlaveDofsType &slaveDofs () const
  {
    slaveDofs_->rebuild();
    return *slaveDofs_;
  }
  class DirichletBuilder;
};

template < class Model, class Space >
class DirichletConstraints< Model,Space > :: DirichletBuilder
    : public CommDataHandleIF< DirichletBuilder, int >
{
public:
  typedef Space SpaceType;
  typedef typename SpaceType::BlockMapperType MapperType;

  enum { nCodim = SpaceType :: GridType :: dimension + 1 };

public:
  typedef int DataType;

  const int myRank_;
  const int mySize_;

  typedef DirichletConstraints< Model,Space > DirichletType;
  const DirichletType &dirichlet_;

  const SpaceType &space_;
  const MapperType &mapper_;

  static const int blockSize = SpaceType::localBlockSize;

public:
  DirichletBuilder( const DirichletType &dirichlet,
                    const SpaceType &space,
                    const MapperType &mapper )
  : myRank_( space.gridPart().comm().rank() ),
    mySize_( space.gridPart().comm().size() ),
    dirichlet_( dirichlet ),
    space_( space ),
    mapper_( mapper )
  {
  }
  bool contains ( int dim, int codim ) const
  {
    return mapper_.contains( codim );
  }

  bool fixedsize ( int dim, int codim ) const
  {
    return false;
  }

  //! read buffer and apply operation
  template< class MessageBuffer, class Entity >
  inline void gather ( MessageBuffer &buffer,
                       const Entity &entity ) const
  {
    unsigned int localBlocks = mapper_.numEntityDofs( entity );
    std::vector<std::size_t> globalBlockDofs(localBlocks);
    mapper_.mapEntityDofs( entity, globalBlockDofs );
    assert( size(entity) == globalBlockDofs.size()*blockSize );
    for( unsigned int localBlock = 0 ; localBlock < globalBlockDofs.size(); ++ localBlock )
    {
      int global = globalBlockDofs[localBlock];
      for (int r=0;r<blockSize;++r)
        if (dirichlet_.dirichletBlocks_[ global ][r] )
          buffer.write( 1 );
        else
          buffer.write( 0 );
    }
  }

  //! read buffer and apply operation
  //! scatter is called for one every entity
  //! several times depending on how much data
  //! was gathered
  template< class MessageBuffer, class EntityType >
  inline void scatter ( MessageBuffer &buffer,
                        const EntityType &entity,
                        size_t n )
  {
    unsigned int localBlocks = mapper_.numEntityDofs( entity );
    std::vector<std::size_t> globalBlockDofs(localBlocks);
    mapper_.mapEntityDofs( entity, globalBlockDofs );
    assert( n == globalBlockDofs.size()*blockSize );
    assert( n == size(entity) );
    for( unsigned int localBlock = 0 ; localBlock < globalBlockDofs.size(); ++ localBlock )
    {
      int global = globalBlockDofs[localBlock];
      for (int r=0;r<blockSize;++r)
      {
        int val;
        buffer.read(val);
        if ( !dirichlet_.dirichletBlocks_[ global ][r] && val == 1)
          dirichlet_.dirichletBlocks_[ global ][r] = true;
      }
    }
  }
  //! return local dof size to be communicated
  template< class Entity >
  size_t size ( const Entity &entity ) const
  {
    return blockSize * mapper_.numEntityDofs( entity );
  }
};

} // end namespace Dune
#endif
