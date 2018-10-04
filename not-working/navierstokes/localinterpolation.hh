#ifndef SPACE_LOCALINTERPOLATION_HH
#define SPACE_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{

  namespace Fem
  {

    template< class FunctionSpace >
    struct LocalBubbleElementInterpolation
    {
      typedef typename FunctionSpace::DomainType DomainType;
      typedef typename FunctionSpace::RangeType RangeType;
      static const int dimDomain = FunctionSpace::dimDomain;
      static const int dimRange = FunctionSpace::dimRange;

      LocalBubbleElementInterpolation ()
        : points_( dimDomain + 2, DomainType( 0.0 ) )
      {
        for( int i = 0; i < dimDomain; ++i )
          points_[ i + 1 ][ i ] = 1.0;

        points_[ dimDomain +1 ] = DomainType( 1.0 / ( dimDomain + 1.0 ) );
      }

      LocalBubbleElementInterpolation ( const LocalBubbleElementInterpolation & ) = default;
      LocalBubbleElementInterpolation ( LocalBubbleElementInterpolation && ) = default;

      //! [Evaluation of local interpolations]
      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &lf, LocalDofVector &ldv ) const
      //! [Evaluation of local interpolations]
      {
        int k = 0;
        for( const DomainType &x : points_ )
        {
          RangeType phi;
          lf.evaluate( x, phi );
          for( int i = 0; i < dimRange; ++i )
            ldv[ k++ ] = phi[ i ];
        }
      }

    private:
      std::vector< DomainType > points_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef SPACE_LOCALINTERPOLATION_HH
