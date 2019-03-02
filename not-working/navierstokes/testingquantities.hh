#ifndef TESTINGQUANTITIES_HH
#define TESTINGQUANTITIES_HH

// iostream includes
#include <iostream>
#include <string>
#include <sstream>

//this function outputs the cahn-hilliard energy of the system created from the full solution.
template <class ExactDiscreteFunction, class DiscreteFunction, class ProblemType>
double* testintegral(ExactDiscreteFunction &exactFinal,DiscreteFunction &uFinal, ProblemType &problem )
{
  typedef typename DiscreteFunction :: DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename ExactDiscreteFunction::LocalFunctionType ExactLocalFunctionType;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
  // typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
  //typedef typename DiscreteFunctionType :: DomainType DomainType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
  // typedef typename GeometryType :: LocalCoordinate LocalCoordinateType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  //typedef Dune::Fem::TimeProviderBase TimeProviderType;

  static double testresults[1];
  const DiscreteFunctionSpaceType &dfSpace(uFinal.space());
  //const double endTime  = Dune::Fem::Parameter::getValue< double >( "navierstokes.endtime", 2.0 );
  //values of interest
  // static const int dimension = GridType :: dimension;
  double lTwoDiff = 0;
 //Loop over the elements
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
    {
      const ElementType &entity = *it;
      // get elements geometry
      const GeometryType &geometry = entity.geometry();
      // get local representation of the discrete functions
      const ExactLocalFunctionType exactLocal = exactFinal.localFunction( entity );
      const LocalFunctionType uLocal = uFinal.localFunction( entity );
      // obtain quadrature order
      const int quadOrder = uLocal.order();
      { // element integral
    QuadratureType quadrature( entity, quadOrder );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {

        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        RangeType exact; RangeType u;
        exactLocal.evaluate( quadrature[ pt ], exact );
        uLocal.evaluate( quadrature[ pt ], u );


        //------------------Velocity L^2 difference
        lTwoDiff += (std::abs(u[0]-exact[0])+std::abs(u[1]-exact[1]))*weight;


        //------FINISH VOLUME INTEGRALS
        // const DomainType xGlobal = entity.geometry().global(coordinate(quadrature[pt]));//This is correct

      }
      }
    }

  testresults[0]=lTwoDiff;
  return testresults;
}
#endif
