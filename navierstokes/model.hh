#ifndef MYDIFFMODEL
#define MYDIFFMODEL

#include <dune/fem/schemes/diffusionmodel.hh>
#include "modelinterface.hh"

template< class FunctionSpace, class GridPart >
struct MyDiffusionModel;
template< class FunctionSpace, class GridPart >
struct MyDiffusionModelTraits
{
  typedef GridPart GridPartType;
  static const int dimRange = FunctionSpace::dimRange;
  typedef MyDiffusionModel<FunctionSpace,GridPart> ModelType;
};


template< class FunctionSpace, class GridPart >
struct MyDiffusionModel : public DiffusionModelEngine<MyDiffusionModelTraits<FunctionSpace,GridPart>>
{
  typedef DiffusionModelEngine<MyDiffusionModelTraits<FunctionSpace,GridPart>> BaseType;
  typedef GridPart GridPartType;
  typedef double RangeFieldType;
  static const int dimRange = FunctionSpace::dimRange;
  static const int dimDomain = GridPartType::dimensionworld;
  static const int dimLocal  = GridPartType::dimension;
  typedef typename GridPart::template Codim<0>::EntityType EntityType;
  typedef typename GridPart::IntersectionType IntersectionType;
  typedef Dune::Fem::FunctionSpace< double, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  typedef ProblemInterface< FunctionSpaceType > ProblemType ;

  MyDiffusionModel( const ProblemType& problem )
    : BaseType(*this), problem_( problem )
  {
  }

  template <class Point>
  void source ( const Point &point,
           const RangeType &value,
           const JacobianRangeType &gradient,
           RangeType &flux ) const
  {
    linSource( value, gradient, point, value, gradient, flux );
  }
  template <class Point>
  void linSource ( const RangeType& uBar,
              const JacobianRangeType &gradientBar,
              const Point &point,
              const RangeType &value,
              const JacobianRangeType &gradient,
              RangeType &flux ) const
  {
    const DomainType xGlobal = entity_->geometry().global( Dune::Fem::coordinate( point ) );
    RangeType m;
    problem_.m(xGlobal,m);
    for (unsigned int i=0;i<flux.size();++i)
      flux[i] = m[i]*value[i];
  }
  template <class Point>
  void diffusiveFlux ( const Point &point,
                  const RangeType &value,
                  const JacobianRangeType &gradient,
                  JacobianRangeType &flux ) const
  {
    linDiffusiveFlux( value, gradient, point, value, gradient, flux );
  }
  template <class Point>
  void linDiffusiveFlux ( const RangeType& uBar,
                     const JacobianRangeType& gradientBar,
                     const Point &point,
                     const RangeType &value,
                     const JacobianRangeType &gradient,
                     JacobianRangeType &flux ) const
  {
    flux = gradient;
  }
  template <class Point>
  void fluxDivergence( const Point &point,
                       const RangeType &value,
                       const JacobianRangeType &jacobian,
                       const HessianRangeType &hessian,
                       RangeType &result) const
  {
    result = RangeType(0);
  }
  template <class Point>
  void alpha(const Point &point,
        const RangeType &value,
        RangeType &val) const
  {
    linAlpha(value,point,value,val);
  }
  template <class Point>
  void linAlpha(const RangeType &uBar,
           const Point &point,
           const RangeType &value,
           RangeType &val) const
  {
    const DomainType xGlobal = entity_->geometry().global( Dune::Fem::coordinate( point ) );
    RangeType alpha;
    problem_.alpha(xGlobal,alpha);
    for (unsigned int i=0;i<val.size();++i)
      val[i] = alpha[i]*value[i];
  }
  bool hasDirichletBoundary () const
  {
    return problem_.hasDirichletBoundary() ;
  }
  bool hasNeumanBoundary () const
  {
    return problem_.hasNeumanBoundary() ;
  }
  bool isDirichletIntersection( const IntersectionType& inter, Dune::FieldVector<int,dimRange> &dirichletComponent ) const
  {
    dirichletComponent = Dune::FieldVector<bool,dimRange>( problem_.isDirichletPoint( inter.geometry().center() ) );
    return dirichletComponent[0];
  }
  template <class Point>
  void dirichlet( int boundaryId,
                  const Point &point,
                  RangeType &value) const
  {
    const DomainType xGlobal = entity_->geometry().global( Dune::Fem::coordinate( point ) );
    problem_.g(xGlobal,value);
  }
  void f(const DomainType& x, RangeType& value) const
  {
    problem_.f(x,value);
  }
  void n(const DomainType& x, RangeType& value) const
  {
    problem_.n(x,value);
  }
  void exact(const DomainType& x, RangeType& value) const
  {
    problem_.g(x,value);
  }
  void jacobianExact(const DomainType& x, JacobianRangeType& value) const
  {
    value = JacobianRangeType(0);
  }
  bool init(const EntityType &entity) const
  {
    entity_ = &entity;
    return true;
  }
  mutable const EntityType *entity_;
  const ProblemType &problem_;
};

#endif // MYDIFFMODEL
