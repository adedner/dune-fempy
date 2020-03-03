#include <iostream>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/function/common/rangegenerators.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/common/bindguard.hh>

// Attention: the coefficients of the functional will be ADDED to RHSStorage!!!
template <class Point, class Functional, class Error>
double pointFunctional(const Point &point, Functional &functional, Error &error)
{
  typedef typename Functional::DiscreteFunctionSpaceType::GridPartType GridPartType;
  Dune::Fem::EntitySearch<GridPartType> search(functional.space().gridPart());
  const auto &entity = search(point);
  const auto localPoint = entity.geometry().local(point);

  Dune::Fem::AddLocalContribution< Functional > wLocal( functional );
  {
    auto guard = Dune::Fem::bindGuard( wLocal, entity );
    typename Functional::RangeType one( 1 );
    wLocal.axpy(localPoint, one);
  }
  auto localError = constLocalFunction(error);
  localError.bind(entity);
  return localError.evaluate(localPoint)[0];
}
