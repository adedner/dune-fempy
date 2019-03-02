#ifndef DUNE_FEM_EDGESPACE_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_EDGESPACE_SHAPEFUNCTIONSET_HH

#include <dune/common/exceptions.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/simple.hh>

namespace Dune
{

  namespace Fem
  {

    // SimplexBubbleElementShapeFunctionSet
    // ----------------------------------

    template< class FunctionSpace >
    class SimplexBubbleElementShapeFunctionSet
    {
      typedef SimplexBubbleElementShapeFunctionSet< FunctionSpace > ThisType;
    public:
      static const int dimDomain =  FunctionSpace :: dimDomain;
      static const int polynomialOrder = dimDomain + 1;
      static const int numShapeFunctions = dimDomain + 2;

      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpace :: DomainType DomainType;
      typedef typename FunctionSpace :: RangeType RangeType;
      static_assert( RangeType::dimension == 1, "This is a scalar shapefunction set" );
      typedef typename FunctionSpace :: JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpace :: HessianRangeType HessianRangeType;

      //! [Main methods for shape functions]
      SimplexBubbleElementShapeFunctionSet () {}

      template<class GeometryType >
      SimplexBubbleElementShapeFunctionSet ( const GeometryType& gt )
      {
        if( !gt.isSimplex() )
          DUNE_THROW( NotImplemented, "BubbleElementShapeFunctionSet is implemented for Simplex only." );
      }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      //! [Main methods for shape functions]
      {
        DomainType xRef = coordinate( x );
        RangeType phi(1), phi0(1);
        for( int i=0; i< dimDomain; ++i )
        {
          functor( i+1, RangeType( xRef[ i ] ) );
          phi0[ 0 ] -= xRef[ i ];
          phi[ 0 ] *= xRef[ i ] ;
        }

        phi[ 0 ] *= phi0[ 0 ] / std::pow( ( dimDomain + 1.0 ), dimDomain + 1.0 );
        functor( 0, phi0 );
        functor( dimDomain +1, phi );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        DomainType xRef = coordinate( x );

        JacobianRangeType jac(0), jac0( -1 );
        RangeType phi0( 1 );

        functor( 0, jac0 );

        for( int i=0; i< dimDomain; ++i )
        {
          phi0[ 0 ] -= xRef[ i ];

          for( int j=1; j < dimDomain; ++j )
            jac0[ 0 ][ (i+j)%dimDomain ] *= xRef[ i ];

          jac[ 0 ][ i ] = 1;
          functor( i+1, jac );
          jac[ 0 ][ i ] = 0;
        }

        for( int i=0; i< dimDomain; ++i )
          jac0[ 0 ][ i ] *= -(phi0[ 0 ] - xRef[ i ]);
        jac0[ 0 ] *= 1.0 / std::pow( dimDomain + 1.0, dimDomain + 1.0 );
        functor( dimDomain +1, jac0 );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        DUNE_THROW( NotImplemented, "NotImplemented" );
        DomainType xRef = coordinate( x );
        HessianRangeType hes;
        functor( 0, hes );
        functor( 1, hes );
        functor( 2, hes );
        functor( 3, hes );
      }

      int order () const { return dimDomain + 1; }

      std::size_t size () const { return dimDomain +2; }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_EDGESPACE_SHAPEFUNCTIONSET_HH
