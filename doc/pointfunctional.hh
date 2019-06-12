#include <vector>
#include <utility>
#include <tuple>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/common/bindguard.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#if 0
template <class Entity, class Error>
double computeError(const Entity &entity, const Error &error)
{
  double ret = 0;
  typedef typename Error::DiscreteFunctionSpaceType::GridPartType GridPartType;
  Dune::Fem::ConstLocalFunction<Error> constLF(entity, error);
  Dune::Fem::CachingQuadrature< GridPartType, 0 > quadrature( entity, 1 );
  const size_t numQuadraturePoints = quadrature.nop();
  for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    ret += constLF.evaluate(quadrature[pt])[0]*quadrature.weight(pt);
  return ret;
}
template <class Point, class Functional, class... Errors>
auto
pointFunctional(const Point &point, Functional &functional, const Errors&... errors)
{
  typedef typename Functional::DiscreteFunctionSpaceType::GridPartType GridPartType;
  Dune::Fem::EntitySearch<GridPartType> search(functional.space().gridPart());
  const auto &entity = search(point);

  functional.clear();
  double weight = 0;
  Dune::Fem::AddLocalContribution< Functional > wLocal( functional );
  for( const auto &entity : functional.space() )
  {
    auto guard = Dune::Fem::bindGuard( wLocal, entity );
    Dune::Fem::CachingQuadrature< GridPartType, 0 > quadrature( entity, 5 );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      std::cout << (point-entity.geometry().global(quadrature.point(pt))).two_norm() << std::endl;
      if ((point-entity.geometry().global(quadrature.point(pt))).two_norm()>0.01) continue;
      typename Functional::RangeType one( quadrature.weight(pt) );
      wLocal.axpy(quadrature[pt], one);
      weight += quadrature.weight(pt);
    }
  }
  for (auto d = functional.dbegin(); d != functional.dend(); ++d)
    *d /= weight;
  return std::make_tuple( computeError(functional.space(),point,errors)... );
}
#else
template <class Entity, class LocalPoint, class Error>
double computeError(const Entity &entity, const LocalPoint &localPoint, const Error &error)
{
  Dune::Fem::ConstLocalFunction<Error> constLF(entity, error);
  return constLF.evaluate(localPoint)[0];
}
#endif

template <class Point, class Functional, class... Errors>
auto
pointFunctional(const Point &point, Functional &functional, const Errors&... errors)
{
  typedef typename Functional::DiscreteFunctionSpaceType::GridPartType GridPartType;
  Dune::Fem::EntitySearch<GridPartType> search(functional.space().gridPart());
  const auto &entity = search(point);
  const auto localPoint = entity.geometry().local(point);

  functional.clear();
  Dune::Fem::AddLocalContribution< Functional > wLocal( functional );
  {
    auto guard = Dune::Fem::bindGuard( wLocal, entity );
    typename Functional::RangeType one( 1 );
    wLocal.axpy(localPoint, one);
  }
  return std::make_tuple( computeError(entity,localPoint,errors)... );
}
