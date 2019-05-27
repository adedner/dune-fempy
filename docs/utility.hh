#include <cstddef>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/fem/space/shapefunctionset/orthonormal.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#if HAVE_DUNE_FEM_DG
#include <dune/fem-dg/operator/limiter/limiter.hh>
#endif

#include <dune/fem/misc/linesegmentsampler.hh>

#if 0
template <class DiscreteFunction>
int minMax( const DiscreteFunction& solution )
{
  typedef typename DiscreteFunction :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType ::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType ::GridPartType GridPartType;
  const DiscreteFunctionSpaceType& space = solution.space();
  const int dimRange = DiscreteFunctionSpaceType :: dimRange;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > Quadrature;

  RangeType minVal( 1e308 );

  bool isNan = false;

  std::vector< RangeType > values;
  for( const auto& element : space )
  {
    Quadrature quad( element, 2*space.order( element )+3 );
    auto lf = solution.localFunction( element );
    const int nop = quad.nop();
    values.resize( nop );
    lf.evaluateQuadrature( quad, values );
    for( int i=0; i<nop; ++i )
    {
      RangeType& val = values[ i ];
      for( int d=0; d<dimRange; ++d )
      {
        minVal[d] = std::min( minVal[d], val[ d ] );
        isNan = !(val[d]==val[d]);
      }
    }
  }
  // std::cout << "Min values: p = " << minVal[0] << "  s = " << minVal[1] << std::endl;
  if( minVal[ 0 ] < 0 || minVal[ 1] < 0 )
    return -1;
  if( isNan )
    return 1;
  return 0;
}
#endif
template <class GF, class DT>
std::pair<std::vector<DT>, std::vector<typename GF::RangeType>>
sample(const GF &gf, DT &start, DT &end, int n)
{
  Dune::Fem::LineSegmentSampler<typename GF::GridPartType> sampler(gf.gridPart(),start,end);
  std::vector<DT> coords(n);
  std::vector<typename GF::RangeType> values(n);
  sampler(gf,values);
  sampler.samplePoints(coords);
  return std::make_pair(coords,values);
}
