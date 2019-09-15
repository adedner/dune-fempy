#include <iostream>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/function/common/rangegenerators.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/common/bindguard.hh>

// Attention: the coefficients of the functional will be ADDED to RHSStorage!!!
template<class RHSStorage, class RHSFunction>
void assembleRHS(RHSStorage& functional, const RHSFunction& function)
{
  const unsigned quadOrder = function.order() + functional.order();

  auto localFunction = constLocalFunction(function);
  Dune::Fem::AddLocalContribution<RHSStorage> localFunctional(functional);
  for (const auto& entity : Dune::Fem::entities(functional)) {
    auto quadrature = elementQuadrature(functional.gridPart(), entity, quadOrder);

    localFunction.bind(entity);
    auto guard = Dune::Fem::bindGuard(localFunctional, entity);

    const auto& geometry = entity.geometry();

    for (std::size_t qp = 0, nop = quadrature.nop(); qp < nop; ++qp) {
      const auto weight = quadrature.weight(qp) * geometry.integrationElement(quadrature.localPoint(qp));
      localFunctional.axpy(quadrature[qp], weight * localFunction.evaluate(quadrature[qp]));
    }
  }
}

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
