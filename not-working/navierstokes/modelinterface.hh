#ifndef ELLIPTC_MODELINTERFACE_HH
#define ELLIPTC_MODELINTERFACE_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>

template< class FunctionSpace, class GridPart >
struct DiffusionModelInterface
{
  typedef FunctionSpace FunctionSpaceType;
  typedef GridPart GridPartType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  static const int dimRange = FunctionSpaceType::dimRange;

  static const bool isLinear = true;
  static const bool isSymmetric = true;

  template< class Entity, class Point >
  void source ( const Entity &entity,
                const Point &x,
                const RangeType &value,
                const JacobianRangeType &gradient,
                RangeType &flux ) const
  {
    linSource( value, gradient, entity, x, value, gradient, flux );
  }

  // the linearization of the source function
  template< class Entity, class Point >
  void linSource ( const RangeType& uBar,
                   const JacobianRangeType &gradientBar,
                   const Entity &entity,
                   const Point &x,
                   const RangeType &value,
                   const JacobianRangeType &gradient,
                   RangeType &flux ) const
  {
    flux = 0;
  }
  //! return the diffusive flux
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity,
                       const Point &x,
                       const RangeType &value,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
    linDiffusiveFlux( value, gradient, entity, x, value, gradient, flux );
  }
  // linearization of diffusiveFlux
  template< class Entity, class Point >
  void linDiffusiveFlux ( const RangeType& uBar,
                          const JacobianRangeType& gradientBar,
                          const Entity &entity,
                          const Point &x,
                          const RangeType &value,
                          const JacobianRangeType &gradient,
                          JacobianRangeType &flux ) const
  {
    // the flux is simply the identity
    flux = 0;
  }
  template< class Entity, class Point >
  void alpha(const Entity &entity, const Point &x,
             const RangeType &value,
             RangeType &val) const
  {
    linAlpha(value,entity,x,value,val);
  }
  template< class Entity, class Point >
  void linAlpha(const RangeType &uBar,
                const Entity &entity, const Point &x,
                const RangeType &value,
                RangeType &val) const
  {
    val = 0;
  }
  //! extract some methods from the problem class for boundary traatment
  bool hasDirichletBoundary () const
  {
    return false;
  }
  bool hasNeumanBoundary () const
  {
    return false;
  }
  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<int,dimRange> &dirichletComponent ) const
  {
    return false;
  }
};

template <class FullDF>
struct LocalVelocityExtractor
{
  typedef typename FullDF::GridPartType GridPartType;
  typedef Dune::Fem::FunctionSpace<double,double, GridPartType::dimensionworld, GridPartType::dimensionworld > FunctionSpaceType;
  typedef Dune::Fem::FunctionSpace<double,double, GridPartType::dimensionworld, GridPartType::dimensionworld+1 > FullFunctionSpaceType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FullDF::EntityType EntityType;

  LocalVelocityExtractor( const FullDF df )
  : df_( df ), ldf_(df) {}

  //! evaluate function
  template <class Point>
  void evaluate( const Point& x, RangeType& ret ) const
  {
    typename FullFunctionSpaceType::RangeType fullRet;
    ldf_.evaluate(x,fullRet);
    for (int i=0;i<RangeType::dimension;++i)
      ret[i] = fullRet[i];
  }
  //! jacobian function (only for exact)
  void jacobian( const DomainType& x, JacobianRangeType& ret ) const
  {
    typename FullFunctionSpaceType::JacobianRangeType fullJacobianRet;
    ldf_.evaluate(x,fullJacobianRet);
    for (int i=0;i<RangeType::dimension;++i)
      ret[i] = fullJacobianRet[i];
  }
  void init(const typename FullDF::EntityType &entity)
  {
    ldf_.init(entity);
  }
  FullDF df_;
  typename FullDF::LocalFunctionType ldf_;
};



#endif // #ifndef ELLIPTC_MODEL_HH
