#ifndef BURGERS_MODEL_HH
#define BURGERS_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/blockvectorfunction.hh>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

//#include "temporalprobleminterface.hh"
#include "navierstokes.hh"
#include "model.hh"

template< class Model >
struct BurgersStateModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                    Model::GridPartType::dimensionworld, Model::GridPartType::dimensionworld > ,typename Model::GridPartType>
{
  typedef typename Model::GridPartType GridPartType;
  typedef Model ModelType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld > VelocityFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, 1 > PressureFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld+1 >FullFunctionSpaceType;
  static const int dimRange = VelocityFunctionSpaceType::dimRange;
  static const int FulldimRange = FullFunctionSpaceType::dimRange;

  typedef typename VelocityFunctionSpaceType::DomainType DomainType;
  typedef typename VelocityFunctionSpaceType::RangeType RangeType;
  typedef typename VelocityFunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename VelocityFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename VelocityFunctionSpaceType::RangeFieldType RangeFieldType;

  //! constructor
  BurgersStateModel( const Model& model,
             const GridPartType &gridPart,
             const double dt,
             const double& alphaOne,
             const double& alphaTwo,
             const bool implicit)
    : model_(model),
      gridPart_(gridPart),
      dt_(dt),
      alphaOne_(alphaOne),
      alphaTwo_(alphaTwo),
      zeroVelocityFunction_( "zero velocity", localZeroVelocity_, gridPart ),
      implicit_(implicit),
      nbc_(model_.neumanBoundary(gridPart)),
      rhs_(model_.rightHandSide(gridPart))
  {
    if (implicit==true)
    {
      timeStepFactor_=1.0;
      timeStepTheta_ = 1.0;}
    else
    {
      timeStepFactor_ = -1.0;
      double theta = Dune::Fem::Parameter::getValue< double >("navierstokes.implicitfactor",0.585786);
      timeStepTheta_ = (1.0-theta)/theta;
    }
  }

  template< class Entity, class Point >
  void source ( const Entity &entity,
                const Point &x,
                const RangeType &value,
                const JacobianRangeType &gradient,
                RangeType &source ) const
  {
    linSource( value, gradient, entity, x, value, gradient, source );
  }

  // the linearization of the source function
  template< class Entity, class Point >
  void linSource ( const RangeType& uBar,
                   const JacobianRangeType& gradientBar,
                   const Entity &entity,
                   const Point &x,
                   const RangeType &value,
                   const JacobianRangeType &gradient,
                   RangeType &source ) const
  {
    source = value;
    source *=  (timeStepFactor_/dt_)*alphaOne_;
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
    flux=gradient;
    flux*= alphaTwo_*timeStepTheta_;

  }

  //! extract some methods from the problem class for boundary traatment
  bool hasDirichletBoundary () const
  {
    return model_.hasDirichletBoundary() ;
  }
  bool hasNeumanBoundary () const
  {
    return model_.hasNeumanBoundary() ;
  }

  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<int,dimRange> &dirichletComponent ) const
  {
    Dune::FieldVector<int,dimRange+1> d;
    bool r = model_.isDirichletIntersection( inter, d );
    for (int i=0;i<dimRange;++i) dirichletComponent[i] = d[i];
    return r;
  }

protected:
  //for dirichlet zero BCS
  typedef GridPartType GP;
  struct LocalZeroVelocityFunction
  {
    typedef GP GridPartType;
    static const int dimDomain = GridPartType::dimensionworld;
    typedef Dune::Fem::FunctionSpace<double,double,dimDomain,dimDomain> FunctionSpaceType;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename VelocityFunctionSpaceType::RangeType RangeType;
    template< class PointType >
    void evaluate ( const PointType &x, RangeType &val )
    {
      val = 0;
    }
    void init ( const EntityType &entity )
    {}
  };
  typedef Dune::Fem::LocalFunctionAdapter< LocalZeroVelocityFunction > ZeroVelocityFunctionType;

public:
  const ZeroVelocityFunctionType &zeroVelocity() const
  {
    return zeroVelocityFunction_;
  }
  Dune::Fem::LocalFunctionAdapter< LocalVelocityExtractor<typename ModelType::NeumanBoundaryType> > neumanBoundary(  ) const
  {
    return Dune::Fem::LocalFunctionAdapter< decltype(nbc_) >("right hand side", nbc_, gridPart_ );
  }
  // return Fem :: Function for right hand side
  Dune::Fem::LocalFunctionAdapter< LocalVelocityExtractor<typename ModelType::RightHandSideType> > rightHandSide(  ) const
  {
    return Dune::Fem::LocalFunctionAdapter< decltype(rhs_) >("right hand side", rhs_, gridPart_ );
  }
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  void init(const EntityType &entity) const
  {
    const_cast<ModelType&>(model_).init(entity);
  }
  template <class Point>
  void dirichlet( int bndId, const Point &x,
                  RangeType &value) const
  {
    typename FullFunctionSpaceType::RangeType fullValue;
    model_.dirichlet(bndId,x,fullValue);
    for (int i=0;i<RangeType::dimension;++i)
      value[i] = fullValue[i];
  }
  class BoundaryWrapper
  {
    const BurgersStateModel<ModelType>& impl_;
    int bndId_;
    public:
    BoundaryWrapper( const BurgersStateModel<ModelType>& impl, int bndId )
    : impl_( impl ), bndId_(bndId) {}

    //! evaluate function
    template <class Point>
    void evaluate( const Point& x, RangeType& ret ) const
    {
      impl_.dirichlet(bndId_,Dune::Fem::coordinate(x),ret);
    }
    //! jacobian function (only for exact)
    void jacobian( const DomainType& x, JacobianRangeType& ret ) const
    {
      DUNE_THROW(Dune::NotImplemented,"rhs jacobian not implemented");
    }
  };

private:
  const Model &model_;
  const GridPartType &gridPart_;
  double dt_;
  const double alphaOne_,alphaTwo_;
  LocalZeroVelocityFunction localZeroVelocity_;
  ZeroVelocityFunctionType zeroVelocityFunction_;
  double timeStepFactor_,timeStepTheta_;
  const bool implicit_;
  mutable LocalVelocityExtractor<typename ModelType::NeumanBoundaryType> nbc_;
  mutable LocalVelocityExtractor<typename ModelType::RightHandSideType> rhs_;
};

//---------------------------------------------------------------//
//---------------------------------------------------------------//
//---------------------------------------------------------------//


template< class Model >
struct BurgersTransportModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                   Model::GridPartType::dimensionworld, Model::GridPartType::dimensionworld > ,typename Model::GridPartType>
{
  typedef typename Model::GridPartType GridPartType;
  typedef Model ModelType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld > VelocityFunctionSpaceType;
  static const int dimRange = VelocityFunctionSpaceType::dimRange;

  typedef typename VelocityFunctionSpaceType::DomainType DomainType;
  typedef typename VelocityFunctionSpaceType::RangeType RangeType;
  typedef typename VelocityFunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename VelocityFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename VelocityFunctionSpaceType::RangeFieldType RangeFieldType;

  typedef Dune::Fem::FunctionSpace< double, double,
             GridPartType::dimensionworld, GridPartType::dimensionworld+1 > FullFunctionSpaceType;

  //! constructor
  BurgersTransportModel( const Model& model,
                   const GridPartType &gridPart)
    : model_(model),
      gridPart_(gridPart)
  {
  }

  template< class Entity, class Point >
  void source ( const Entity &entity,
                const Point &x,
                const RangeType &value,
                const JacobianRangeType &gradient,
                RangeType &source ) const //value = u^n, gradient = \nabla u^n
  {
     gradient.mv(value,source);//grad^T *value = source
  }

  //! return the diffusive flux
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity,
                       const Point &x,
                       const RangeType &value,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
    flux = JacobianRangeType(0);
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
    Dune::FieldVector<int,dimRange+1> d;
    bool r = model_.isDirichletIntersection( inter, d );
    for (int i=0;i<dimRange;++i) dirichletComponent[i] = d[i];
    return r;
  }
  // return Fem :: Function for Dirichlet boundary values

private:
  const Model &model_;
  const GridPartType &gridPart_;

};

//---------------------------------------------------------------//
//---------------------------------------------------------------//
//---------------------------------------------------------------//

template< class Model, class VeloDF >
struct BurgersDescentModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                            Model::GridPartType::dimensionworld, Model::GridPartType::dimensionworld > ,typename Model::GridPartType>
{
  typedef typename Model::GridPartType GridPartType;
  typedef Model ModelType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld > VelocityFunctionSpaceType;
  static const int dimRange = VelocityFunctionSpaceType::dimRange;


  typedef typename VelocityFunctionSpaceType::DomainType DomainType;
  typedef typename VelocityFunctionSpaceType::RangeType RangeType;
  typedef typename VelocityFunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename VelocityFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename VelocityFunctionSpaceType::RangeFieldType RangeFieldType;

  //for xi_
  typedef VeloDF VelocityDiscreteFunctionType;
  typedef typename VelocityDiscreteFunctionType :: LocalFunctionType LocalVelocityFunctionType;
  typedef typename VeloDF::DiscreteFunctionSpaceType VelocityDiscreteFunctionSpaceType;
  //


  typedef Dune::Fem::FunctionSpace< double, double,
             GridPartType::dimensionworld, GridPartType::dimensionworld+1 > FullFunctionSpaceType;
  typedef NavierStokesProblemInterface< FullFunctionSpaceType> ProblemType ;
  typedef Dune::Fem::TimeProviderBase TimeProviderType;
  //! constructor
  BurgersDescentModel( const Model& model,
               const GridPartType& gridPart,
               const VelocityDiscreteFunctionType& xi)
    : model_(model),
      gridPart_(gridPart),
      xi_(xi)
  {
  }

  template< class Entity, class Point >
  void source ( const Entity &entity,
                const Point &x,
                const RangeType &value,
                const JacobianRangeType &gradient,
                RangeType &source ) const //value = u^n, gradient = \nabla u^n
  {
    const LocalVelocityFunctionType xiLocal = xi_.localFunction(entity);
    RangeType xi;
    xiLocal.evaluate(x,xi);
    gradient.mv(xi,source);
  }

  //! return the diffusive flux
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity,
                       const Point &x,
                       const RangeType &value,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
   flux = JacobianRangeType(0);
    const LocalVelocityFunctionType xiLocal = xi_.localFunction(entity);
    RangeType xi;
    xiLocal.evaluate(x,xi);
    static const int dimDomain = RangeType::dimension;
    //create ((w tenosr xi ) : nabla phi
    for( unsigned int localRow = 0; localRow < dimDomain; ++localRow )
      {
    for( unsigned int localCol = 0; localCol < dimDomain; ++localCol )
      {
        flux[localRow][localCol]+= value[localRow]*xi[localCol];//sum_j (sum_i w_i xi_j) Gphi_ij
      }
      }
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
    Dune::FieldVector<int,dimRange+1> d;
    bool r = model_.isDirichletIntersection( inter, d );
    for (int i=0;i<dimRange;++i) dirichletComponent[i] = d[i];
    return r;
  }
  // return Fem :: Function for Dirichlet boundary values

protected:

private:
  const Model &model_;
  const GridPartType &gridPart_;
  const VelocityDiscreteFunctionType &xi_;

};

//---------------------------------------------------------------//
//---------------------------------------------------------------//
//---------------------------------------------------------------//


template< class Model >
struct BurgersGradModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                    Model::GridPartType::dimensionworld, Model::GridPartType::dimensionworld > ,typename Model::GridPartType>
{
  typedef typename Model::GridPartType GridPartType;
  typedef Model ModelType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld > VelocityFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, 1 > PressureFunctionSpaceType;

  typedef typename PressureFunctionSpaceType::DomainType DomainType;

  typedef typename PressureFunctionSpaceType::RangeType DomainRangeType;
  typedef typename PressureFunctionSpaceType::JacobianRangeType DomainJacobianRangeType;
  typedef typename VelocityFunctionSpaceType::RangeType  RangeRangeType;
  typedef typename VelocityFunctionSpaceType::JacobianRangeType RangeJacobianRangeType;

  static const int dimRange = VelocityFunctionSpaceType::dimRange;

  //! constructor
  BurgersGradModel(  const GridPartType &gridPart )
  {
  }

  template< class Entity, class Point >
  void source ( const Entity &entity,
                const Point &x,
                const DomainRangeType &value,
                const DomainJacobianRangeType &gradient,
                RangeRangeType &source ) const
  {
    linSource( value, gradient, entity, x, value, gradient, source );
  }

  // the linearization of the source function
  template< class Entity, class Point >
  void linSource ( const DomainRangeType& uBar,
                   const DomainJacobianRangeType& gradientBar,
                   const Entity &entity,
                   const Point &x,
                   const DomainRangeType &value,
                   const DomainJacobianRangeType &gradient,
                   RangeRangeType &source ) const
  {
    source = 0;
  }

  //! return the diffusive flux
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity,
                       const Point &x,
                       const DomainRangeType &value,
                       const DomainJacobianRangeType &gradient,
                       RangeJacobianRangeType &flux ) const
  {
    linDiffusiveFlux( value, gradient, entity, x, value, gradient, flux );
  }

  // linearization of diffusiveFlux
  template< class Entity, class Point >
  void linDiffusiveFlux ( const DomainRangeType& uBar,
                          const DomainJacobianRangeType& gradientBar,
                          const Entity &entity,
                          const Point &x,
                          const DomainRangeType &value,
                          const DomainJacobianRangeType &gradient,
                          RangeJacobianRangeType &flux ) const
  {
    flux = 0.0;

    static const int dimDomain = DomainType::dimension;

    for( unsigned int localRow = 0; localRow < dimDomain; ++localRow )
      flux[ localRow ][ localRow ] -= value[ 0 ]; // -p Id
  }
 template< class Entity, class Point >
  void alpha(const Entity &entity, const Point &x,
             const DomainRangeType &value,
             RangeRangeType &val) const
  {
    linAlpha(value,entity,x,value,val);
  }
  template< class Entity, class Point >
  void linAlpha(const DomainRangeType &uBar,
                const Entity &entity, const Point &x,
                const DomainRangeType &value,
                RangeRangeType &val) const
  {
    val = 0;
  }
  bool hasDirichletBoundary () const
  {
    return false;
  }
  bool hasNeumanBoundary () const
  {
    return false;
  }
};
#endif // #ifndef BURGERS_MODEL_HH
