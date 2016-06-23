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

template< class GridPart >
struct BurgersStateModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                    GridPart::dimensionworld, GridPart::dimensionworld > ,GridPart>
{
  typedef GridPart GridPartType;
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

  typedef NavierStokesProblemInterface< FullFunctionSpaceType > ProblemType ;
  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  //BCs
protected:
  enum FunctionId { rhs, bndD, bndN };
  template <FunctionId id>
  class FunctionWrapper;
public:
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<rhs>, GridPartType > RightHandSideType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<bndD>, GridPartType > DirichletBoundaryType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<bndN>, GridPartType > NeumanBoundaryType;



  //! constructor
  BurgersStateModel( const ProblemType& problem,
             const GridPartType &gridPart,
             const double& alphaOne,
             const double& alphaTwo,
             const bool implicit)
    : problem_(problem),
      timeProvider_(problem.timeProvider()),
      gridPart_(gridPart),
      zeroVelocityFunction_( "zero velocity", localZeroVelocity_, gridPart ),
      rhs_(problem_),
      bndD_(problem_),
      bndN_(problem_),
      alphaOne_(alphaOne),
      alphaTwo_(alphaTwo),
      implicit_(implicit)//boolean
  {
    if (implicit==true){
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
    source *=  (timeStepFactor_/timeProvider_.deltaT())*alphaOne_;
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
    return problem_.hasDirichletBoundary() ;
  }
  bool hasNeumanBoundary () const
  {
    return problem_.hasNeumanBoundary() ;
  }

  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<bool,dimRange> &dirichletComponent ) const
  {
    return problem_.isDirichletPoint( inter.geometry().center() );
  }
  // return Fem :: Function for Dirichlet boundary values
  DirichletBoundaryType dirichletBoundary( ) const
  {
    return DirichletBoundaryType( "boundary function", bndD_, gridPart_, 5 );
  }
  NeumanBoundaryType neumanBoundary( ) const
  {
    return NeumanBoundaryType( "boundary function", bndN_, gridPart_, 5 );
  }
  // return Fem :: Function for right hand side
  RightHandSideType rightHandSide(  ) const
  {
    return RightHandSideType( "right hand side", rhs_, gridPart_, 5 );
  }

protected:
  //for dirichlet zero BCS
 struct LocalZeroVelocityFunction
  {
    typedef GridPart GridPartType;
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

  template <FunctionId id>
  class FunctionWrapper : public Dune::Fem::Function<VelocityFunctionSpaceType, FunctionWrapper< id > >
  {
    const NavierStokesProblemInterface<FullFunctionSpaceType>& impl_;
  public:
    FunctionWrapper( const NavierStokesProblemInterface<FullFunctionSpaceType>& impl )
      : impl_( impl ) {}

    //! evaluate function
    void evaluate( const DomainType& x, RangeType& ret ) const
    {
      typename FullFunctionSpaceType::RangeType fullRet;
      if( id == rhs )
      {
        // call right hand side of implementation
        impl_.f( x, fullRet );
      }
      else if( id == bndD )
      {
        // call dirichlet boudary data of implementation
        impl_.g( x, fullRet );
      }
      else if( id == bndN )
      {
        // call neumann boudary data of implementation
        impl_.n( x, fullRet );
      }
      else
      {
        DUNE_THROW(Dune::NotImplemented,"FunctionId not implemented");
      }
      for (unsigned int i=0;i<ret.size();++i)
        ret[i] = fullRet[i];
    }
  };
public:
  const ZeroVelocityFunctionType &zeroVelocity() const
  {
    return zeroVelocityFunction_;
  }

private:
  const ProblemType &problem_;
  const TimeProviderType &timeProvider_;
  const GridPart &gridPart_;
  LocalZeroVelocityFunction localZeroVelocity_;
  ZeroVelocityFunctionType zeroVelocityFunction_;
  FunctionWrapper<rhs> rhs_;
  FunctionWrapper<bndD> bndD_;
  FunctionWrapper<bndN> bndN_;
  const double alphaOne_,alphaTwo_;
  double timeStepFactor_,timeStepTheta_;
  const bool implicit_;

};

//---------------------------------------------------------------//
//---------------------------------------------------------------//
//---------------------------------------------------------------//


template< class GridPart >
struct BurgersTransportModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                   GridPart::dimensionworld, GridPart::dimensionworld > ,GridPart>
{
  typedef GridPart GridPartType;
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
  typedef NavierStokesProblemInterface< FullFunctionSpaceType > ProblemType ;
  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  //! constructor
  BurgersTransportModel( const ProblemType& problem,
                   const GridPartType &gridPart)
    : problem_(problem),
      timeProvider_(problem.timeProvider()),
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
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<bool,dimRange> &dirichletComponent ) const
  {
    return problem_.isDirichletPoint( inter.geometry().center() );
  }
  // return Fem :: Function for Dirichlet boundary values

protected:

private:
  const ProblemType &problem_;
  const TimeProviderType &timeProvider_;
  const GridPart &gridPart_;

};

//---------------------------------------------------------------//
//---------------------------------------------------------------//
//---------------------------------------------------------------//

template< class GridPart >
struct BurgersDescentModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                    GridPart::dimensionworld, GridPart::dimensionworld > ,GridPart>
{
  typedef GridPart GridPartType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld > VelocityFunctionSpaceType;
  static const int dimRange = VelocityFunctionSpaceType::dimRange;


  typedef typename VelocityFunctionSpaceType::DomainType DomainType;
  typedef typename VelocityFunctionSpaceType::RangeType RangeType;
  typedef typename VelocityFunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename VelocityFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename VelocityFunctionSpaceType::RangeFieldType RangeFieldType;

  //for xi_
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< VelocityFunctionSpaceType, GridPartType, POLORDER+1 > VelocityDiscreteFunctionSpaceType;
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< VelocityDiscreteFunctionSpaceType > VelocityDiscreteFunctionType;
  typedef typename VelocityDiscreteFunctionType :: LocalFunctionType LocalVelocityFunctionType;
  //


  typedef Dune::Fem::FunctionSpace< double, double,
             GridPartType::dimensionworld, GridPartType::dimensionworld+1 > FullFunctionSpaceType;
  typedef NavierStokesProblemInterface< FullFunctionSpaceType> ProblemType ;
  typedef Dune::Fem::TimeProviderBase TimeProviderType;
  //! constructor
  BurgersDescentModel( const ProblemType& problem,
               const GridPartType& gridPart,
               const VelocityDiscreteFunctionType& xi)
    : problem_(problem),
      timeProvider_(problem.timeProvider()),
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
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<bool,dimRange> &dirichletComponent ) const
  {
    return problem_.isDirichletPoint( inter.geometry().center() );
  }
  // return Fem :: Function for Dirichlet boundary values

protected:

private:
  const ProblemType &problem_;
  const TimeProviderType &timeProvider_;
  const GridPart &gridPart_;
  const VelocityDiscreteFunctionType &xi_;

};

//---------------------------------------------------------------//
//---------------------------------------------------------------//
//---------------------------------------------------------------//


template< class GridPart >
struct BurgersGradModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                    GridPart::dimensionworld, GridPart::dimensionworld > ,GridPart>
{
  typedef GridPart GridPartType;
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
