#ifndef DUNE_NOSLIPCONSTRAINTS_HH
#define DUNE_NOSLIPCONSTRAINTS_HH

#include <dune/fem/space/combinedspace/combineddiscretefunctionspace.hh>

#include <dune/fem/function/localfunction/converter.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

#include "dirichletconstraints.hh"

namespace Dune {

  /**A class modelling no-slip constraints for flow problems. It is
   * implicitly assumed that the DiscreteFunctionSpace combines
   * velocity and pressure such that the first Space::dimRange-1
   * dimensions "belong" to the belocity and the last dimension is for
   * the pressure. We provide one version for equal order
   * discretizations for Lagrange-Spaces where the local block size
   * coincides with dimRange and a specialization for the
   * CombinedDiscreteFunctionSpace where the local block size has to
   * be 1.
   */
  template < class Model, class DiscreteFunctionSpace >
  class NoSlipConstraints
    : public DirichletConstraints<Model, DiscreteFunctionSpace>
  {
    typedef DirichletConstraints<Model, DiscreteFunctionSpace> BaseType;
   public:
    typedef Model ModelType;
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

    NoSlipConstraints( const ModelType &model, const DiscreteFunctionSpaceType& space )
      : BaseType(model, space)
    {}
  };

  /**Special implementation for the CombinedDiscreteFunctionSpace. We
   * rather clone the source code from DirichletContraints as somehow
   * DirichletConstraints needs block-size == dimRange all over the
   * place.
   */
  template<class Model, class DFunctionSpace1, class DFunctionSpace2>
  class NoSlipConstraints<Model, Fem::TupleDiscreteFunctionSpace<DFunctionSpace1, DFunctionSpace2> >
    : public DirichletConstraints<Model, DFunctionSpace1>
  {
    typedef DirichletConstraints<Model, DFunctionSpace1> BaseType;
   protected:
    using BaseType::updateDirichletDofs;
    using BaseType::hasDirichletDofs_;
    using BaseType::dirichletDofTreatment;
    using BaseType::space_;

    typedef typename Fem::TupleDiscreteFunctionSpace< DFunctionSpace1, DFunctionSpace2 >::BasisFunctionSetType::template SubBasisFunctionSet< 0 >::type
      SubBasisFunctionSetType;

    struct VelocityConverter
    {
      template< class T, int i >
      FieldVector< T, i-1 > operator() ( const FieldVector< T, i> & in ) const
      {
        FieldVector< T, i-1 > out;
        convert( in, out );
        return out;
      }

      template< class T, int i, int j >
      FieldMatrix< T, i-1, j > operator() ( const FieldMatrix< T, i, j > & in ) const
      {
        FieldVector< T, i-1 > out;
        convert( in, out );
        return out;
      }

    private:
      template< class In , class Out >
      void convert ( const In &in, Out &out ) const
      {
        for( unsigned int i = 0; i< out.size(); ++i )
          out[ i ] = in[ i ];
      }
    };

   public:
    typedef Model ModelType;
    typedef Fem::TupleDiscreteFunctionSpace<DFunctionSpace1, DFunctionSpace2> DiscreteFunctionSpaceType;

    NoSlipConstraints( const ModelType &model, const DiscreteFunctionSpaceType& space )
      : BaseType( model, std::get<0>( space.spaceTuple() ) )
    {
    }

    template< class GridFunctionType, class DiscreteFunctionType >
    void operator() ( const GridFunctionType &u, DiscreteFunctionType &w ) const
    {
      updateDirichletDofs();

      if( hasDirichletDofs_ )
      {
        typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
        typedef typename IteratorType :: Entity EntityType;

        typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
        typedef typename LocalFunctionType::LocalDofVectorType LocalDofVectorType;

        typedef Fem::DenseSubVector< LocalDofVectorType > LocalSubDofVectorType;
        typedef Fem::LocalFunction< SubBasisFunctionSetType, LocalSubDofVectorType > LocalSubFunctionType;

        for( const EntityType &entity : space_ )
        {
          // we still need a reference to the object
          const auto uSubLocal = Fem::localFunctionConverter( u.localFunction( entity ), VelocityConverter() );

          LocalFunctionType wLocal = w.localFunction( entity );
          const SubBasisFunctionSetType subBasisSet = w.space().basisFunctionSet( entity ).template subBasisFunctionSet<0>();
          LocalSubDofVectorType subDofVector( wLocal.localDofVector(), subBasisSet.size(), 0 );
          LocalSubFunctionType wSubLocal( space_.basisFunctionSet( entity ), subDofVector );

          BaseType::dirichletDofTreatment( uSubLocal, wSubLocal );
        }
      }
    }
  };

} // end namespace Dune

#endif
