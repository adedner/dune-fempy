#ifndef Stokes_MODEL_HH
#define Stokes_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

//#include "temporalprobleminterface.hh"
#include "navierstokes.hh"
#include "model.hh"

template< class FunctionSpace, class GridPart >
struct StokesModel : public DiffusionModel<FunctionSpace,GridPart>
{
  typedef DiffusionModel<FunctionSpace,GridPart> BaseType;
  typedef FunctionSpace FunctionSpaceType;
  typedef GridPart GridPartType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld+1 > FullFunctionSpaceType;
typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld > VelocityFunctionSpaceType;
  typedef NavierStokesProblemInterface< FullFunctionSpaceType > ProblemType ;
  typedef Dune::Fem::TimeProviderBase TimeProviderType;
  static const int dimRange = FunctionSpaceType::dimRange;

  //! constructor
  StokesModel( const ProblemType& problem,
             const GridPart &gridPart,
             double mu, double nu)
    : BaseType(problem,gridPart),
      stab_( Dune::Fem::Parameter::getValue< double >( "stokes.stability", 0.0 ) ),
      mu_(mu), nu_(nu)
  {
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
    source = 0;
    for (unsigned int i=0;i<GridPart::dimensionworld;++i)
      source[i] = nu_*value[i];

    static const int dimDomain = RangeType::dimension-1;
    source[dimDomain] = 0;
    for( unsigned int localRow = 0; localRow < dimDomain; ++localRow )
    {
      // conventional incompressibility condition =  nabla dot u times q in dimdomain(=pressure)  row of aphi
      source[ dimDomain ] -= gradient[ localRow ][ localRow ];
    }
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
    flux = 0.0;

    static const int dimDomain = RangeType::dimension-1;

    const auto &geometry = entity.geometry();
    const double scaling = std::pow( geometry.volume() , 1.0 / 3.0 );

    for( unsigned int localRow = 0; localRow < dimDomain; ++localRow )
    {
      // velocity equations -- localRow = 0, 1, 2 = dimDomain - 1
      // symmetric tensor
      for( unsigned int localCol = 0; localCol < dimDomain; ++localCol )
        flux[ localRow ][localCol] += mu_*(gradient[ localRow ][ localCol ] + gradient[ localCol ][ localRow ]);
      // laplace
      // flux[localRow] += mu_*gradient[localRow];
      flux[ localRow ][ localRow ] -= value[ dimDomain ]; // -p Id
      flux[ localRow ][ dimDomain ] = 0;

      // pressure equation -- localRow = dimdomain

      // stabilization term = - k * nabla p dot nabla q  where p and q are pressure shape functions
      flux[ dimDomain ][ localRow ] -= stab_ * scaling * gradient[ dimDomain ][ localRow ];
    }
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<bool,dimRange> &dirichletComponent ) const
  {
    bool isDirichlet = BaseType::isDirichletIntersection(inter, dirichletComponent);
    dirichletComponent[dimRange - 1] = false;

    return isDirichlet;
  }
protected:
  using BaseType::problem_;
  double stab_;
  double mu_,nu_;

};

//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------

template< class GridPart >
struct StokesMainModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                    GridPart::dimensionworld, GridPart::dimensionworld > ,GridPart>
{
  typedef GridPart GridPartType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld > FunctionSpaceType;
  static const int dimRange = FunctionSpaceType::dimRange;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, 1 > PressureFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld+1 > FullFunctionSpaceType;
  typedef NavierStokesProblemInterface< FullFunctionSpaceType > ProblemType ;
  typedef Dune::Fem::TimeProviderBase TimeProviderType;
protected:
  enum FunctionId { rhs, rhsnext, bndD, bndN };
  template <FunctionId id>
  class FunctionWrapper;
public:
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<rhs>, GridPartType > RightHandSideType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<rhsnext>, GridPartType > RightHandSideNextType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<bndD>, GridPartType > DirichletBoundaryType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<bndN>, GridPartType > NeumanBoundaryType;

  //! constructor
  StokesMainModel( const ProblemType& problem,
                   const GridPartType &gridPart,
                   double mu, double nu,const bool implicit)
    : problem_(problem),
      timeProvider_(problem.timeProvider()),
      gridPart_(gridPart),
      localInitialPressure_(problem),
      localInitialVelocity_(problem),
      initialPressure_("initial pressure", localInitialPressure_,gridPart),
      initialVelocity_("initial velocity",localInitialVelocity_,gridPart),
      zeroVelocityFunction_( "zero velocity", localZeroVelocity_, gridPart ),
      rhs_(problem_),
      rhsnext_(problem_),
      bndD_(problem_),
      bndN_(problem_),
      mu_(mu), nu_(nu),
      implicit_(implicit)
  {

    if (implicit==true){
      timeStepFactor_=1.0;
      timeStepTheta_ =1.0;}
    else
      {
    double theta = Dune::Fem::Parameter::getValue< double >("navierstokes.implicitfactor",0.585786);
    timeStepFactor_ = -1.0;
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
    source = 0;
    for (unsigned int i=0;i<GridPart::dimensionworld;++i)
      source[i] = (1./timeProvider_.deltaT())*nu_*value[i];
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
    flux = 0.0;

    static const int dimDomain = RangeType::dimension;

    for( unsigned int localRow = 0; localRow < dimDomain; ++localRow )
    {
      // velocity equations -- localRow = 0, 1, 2 = dimDomain - 1
      // symmetric tensor
      for( unsigned int localCol = 0; localCol < dimDomain; ++localCol )
        flux[ localRow ][localCol] += timeStepFactor_*timeStepTheta_*mu_*(gradient[ localRow ][ localCol ] + gradient[ localCol ][ localRow ]);
      // laplace
      // flux[localRow] +=  timeStepFactor_*timeStepTheta_*mu_*gradient[localRow];
    }
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
    val = RangeType(0);
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


  //-------------------------------------------------------------------------------------
protected:
  class LocalInitialPressureFunction
  {
  public:
    typedef GridPart GridPartType;
    static const int dimDomain = GridPartType::dimensionworld;
    typedef Dune::Fem::FunctionSpace<double,double,dimDomain,dimDomain+1> FullFunctionSpaceType;
    typedef typename FullFunctionSpaceType::RangeType FullRangeType;
    typedef Dune::Fem::FunctionSpace<double,double,dimDomain,dimDomain> VelocityFunctionSpaceType;
    typedef Dune::Fem::FunctionSpace<double,double,dimDomain,1> FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
     typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

    // constructor
    LocalInitialPressureFunction ( const NavierStokesProblemInterface<FullFunctionSpaceType>& problem)
      :    problem_(problem){}

    // evaluate local function
     template< class PointType >
    void evaluate ( const PointType &x, RangeType &val )
    {
      const DomainType xDomain(coordinate(x));
      FullRangeType retP;
      problem_.u(entity_.geometry().global(xDomain),retP);
      val = retP[retP.size()-1];
    }
    void init ( const EntityType &entity )
    {
      entity_=entity;
    }
  private:
    EntityType entity_;
    const NavierStokesProblemInterface<FullFunctionSpaceType>& problem_;
  };

  class LocalInitialVelocityFunction
  {
  public:
    typedef GridPart GridPartType;
    static const int dimDomain = GridPartType::dimensionworld;
    typedef Dune::Fem::FunctionSpace<double,double,dimDomain,dimDomain+1> FullFunctionSpaceType;
    typedef typename FullFunctionSpaceType::RangeType FullRangeType;
    typedef Dune::Fem::FunctionSpace<double,double,dimDomain,dimDomain> FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

    // constructor
    LocalInitialVelocityFunction ( const NavierStokesProblemInterface<FullFunctionSpaceType>& problem)
      :    problem_(problem){}

    // evaluate local function
    template< class PointType >
    void evaluate ( const PointType &x, RangeType &val )
    {
      const DomainType xDomain(coordinate(x));
      FullRangeType retU;
      problem_.u(entity_.geometry().global(xDomain),retU);
      for (unsigned int i=0;i<(retU.size()-1);++i)
    val[i] = retU[i];
    }
    void init ( const EntityType &entity )
    {
      entity_=entity;
    }
  private:
    EntityType entity_;
    const NavierStokesProblemInterface<FullFunctionSpaceType>& problem_;
  };

  struct LocalZeroVelocityFunction
  {
    typedef GridPart GridPartType;
    static const int dimDomain = GridPartType::dimensionworld;
    typedef Dune::Fem::FunctionSpace<double,double,dimDomain,dimDomain> FunctionSpaceType;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename FunctionSpaceType::RangeType RangeType;
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
  class FunctionWrapper : public Dune::Fem::Function< FunctionSpaceType, FunctionWrapper< id > >
  {
    typedef Dune::Fem::TimeProviderBase TimeProviderType;
    public:
    FunctionWrapper( const NavierStokesProblemInterface<FullFunctionSpaceType>& impl)
      : impl_( impl ){}

    //! evaluate function
    void evaluate( const DomainType& x, RangeType& ret ) const
    {
      typename FullFunctionSpaceType::RangeType fullRet;
      if( id == rhs )
      {
        // call right hand side of implementation
        impl_.f( x, fullRet );
      }
      else if( id == rhsnext )
      {
        impl_.fnext( x, fullRet );
      }
      else if( id == bndD )
      {
        // call dirichlet boudary data of implementation
        impl_.g( x, fullRet );
      }
      else if( id == bndN )
      {
        // call dirichlet boudary data of implementation
        impl_.n( x, fullRet );
      }
      else
      {
        DUNE_THROW(Dune::NotImplemented,"FunctionId not implemented");
      }
      for (unsigned int i=0;i<ret.size();++i)
        ret[i] = fullRet[i];
    }
  private:
    const NavierStokesProblemInterface<FullFunctionSpaceType>& impl_;
  };

public:
  typedef Dune::Fem::LocalFunctionAdapter< LocalInitialPressureFunction > InitialPressureFunctionType;
  typedef Dune::Fem::LocalFunctionAdapter< LocalInitialVelocityFunction > InitialVelocityFunctionType;

  const InitialPressureFunctionType &initialPressureFunction() const
  {
    return initialPressure_;
  }
 const InitialVelocityFunctionType &initialVelocityFunction() const
  {
    return initialVelocity_;
  }
  const ZeroVelocityFunctionType &zeroVelocity() const
  {
    return zeroVelocityFunction_;
  }
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
  RightHandSideNextType rightHandSideNext(  ) const
  {
    return RightHandSideNextType( "right hand side - next timestep", rhsnext_, gridPart_, 5 );
  }
  //------------------------------------------------------------------------------------------

private:
  const ProblemType &problem_;
  const TimeProviderType &timeProvider_;
  const GridPart &gridPart_;
  LocalInitialPressureFunction localInitialPressure_;
  LocalInitialVelocityFunction localInitialVelocity_;
  InitialPressureFunctionType initialPressure_;
  InitialVelocityFunctionType initialVelocity_;
  LocalZeroVelocityFunction localZeroVelocity_;
  ZeroVelocityFunctionType zeroVelocityFunction_;
  FunctionWrapper<rhs> rhs_;
  FunctionWrapper<rhsnext> rhsnext_;
  FunctionWrapper<bndD> bndD_;
  FunctionWrapper<bndN> bndN_;
  double mu_,nu_;
  double timeStepFactor_,timeStepTheta_;
  const bool implicit_;
};

//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------

template< class GridPart >
struct StokesGradModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
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
  StokesGradModel(  const GridPartType &gridPart )
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
template< class GridPart >
struct StokesMassModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                    GridPart::dimensionworld, 1 > ,GridPart>
{
  typedef GridPart GridPartType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, 1 > PressureFunctionSpaceType;

  typedef typename PressureFunctionSpaceType::DomainType DomainType;

  typedef typename PressureFunctionSpaceType::RangeType DomainRangeType;
  typedef typename PressureFunctionSpaceType::JacobianRangeType DomainJacobianRangeType;
  typedef typename PressureFunctionSpaceType::RangeType  RangeRangeType;
  typedef typename PressureFunctionSpaceType::JacobianRangeType RangeJacobianRangeType;

  static const int dimRange = PressureFunctionSpaceType::dimRange;

  //! constructor
  StokesMassModel(  const GridPartType &gridPart )
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
    source = value;
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
template< class GridPart >
struct StokesDivergenceModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                    GridPart::dimensionworld, 1 > ,GridPart>
{
  typedef GridPart GridPartType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld > VelocityFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, 1 > PressureFunctionSpaceType;

  typedef typename VelocityFunctionSpaceType::DomainType DomainType;

  typedef typename VelocityFunctionSpaceType::RangeType  DomainRangeType;
  typedef typename VelocityFunctionSpaceType::JacobianRangeType DomainJacobianRangeType;
  typedef typename PressureFunctionSpaceType::RangeType RangeRangeType;
  typedef typename PressureFunctionSpaceType::JacobianRangeType RangeJacobianRangeType;

  static const int dimRange = PressureFunctionSpaceType::dimRange;

  //! constructor
  StokesDivergenceModel(  const GridPartType &gridPart )
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
    static const int dimDomain = DomainType::dimension;
    for( unsigned int localRow = 0; localRow < dimDomain; ++localRow )
      source[ 0 ] += gradient[ localRow ][ localRow ];
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
    flux = 0;
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
template< class GridPart >
struct StokesPrecondModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,GridPart::dimensionworld,1> ,GridPart>
{
  typedef GridPart GridPartType;
  typedef Dune::Fem::FunctionSpace< double, double,GridPartType::dimensionworld,1> FunctionSpaceType;
  static const int dimRange = FunctionSpaceType::dimRange;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld+1 > FullFunctionSpaceType;
 typedef Dune::Fem::FunctionSpace< double, double,
              GridPartType::dimensionworld, GridPartType::dimensionworld > VelocityFunctionSpaceType;
   typedef NavierStokesProblemInterface< FullFunctionSpaceType> ProblemType ;
  typedef Dune::Fem::TimeProviderBase TimeProviderType;

public:
  //! constructor
  StokesPrecondModel( const ProblemType& problem,
                      const GridPartType &gridPart,
                      double mu, double nu)
    : problem_(problem),
      gridPart_(gridPart),
      zeroVelocityFunction_( "zero velocity", localZeroVelocity_, gridPart ),
      mu_(mu),
      nu_(nu)
  {
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
    source = 0;
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
    flux = gradient;
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
    // robin b.c on boundaries that are not dirichlet
    if ( !problem_.isDirichletPoint( entity.geometry().global(coordinate(x)) ) )
      for (unsigned int i=0;i<val.size();++i)
        val[i] = nu_*value[i];
  }
  //! extract some methods from the problem class for boundary traatment
  bool hasDirichletBoundary () const
  {
    return false;
  }
  bool hasNeumanBoundary () const
  {
    return true;
  }

  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<bool,dimRange> &dirichletComponent ) const
  {
    return false;
  }
protected:
  struct LocalZeroVelocityFunction
  {
    typedef GridPart GridPartType;
    static const int dimDomain = GridPartType::dimensionworld;
    typedef Dune::Fem::FunctionSpace<double,double,dimDomain,1> FunctionSpaceType;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename FunctionSpaceType::RangeType RangeType;
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
private:
  const ProblemType &problem_;
  const GridPart &gridPart_;
  LocalZeroVelocityFunction localZeroVelocity_;
  ZeroVelocityFunctionType zeroVelocityFunction_;
  double mu_,nu_;
};

//---------------------------------------------------------------//
//---------------------------------------------------------------//
//---------------------------------------------------------------//


template< class GridPart >
struct StokesTransportModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
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
  typedef ProblemType InitialFunctionType;
  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  //! constructor
  StokesTransportModel(const ProblemType &problem, const GridPartType &gridPart)
  : gridPart_(gridPart), problem_(problem)
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
  const InitialFunctionType &initialFunction() const
  {
    return problem_;
  }
protected:

private:
  const GridPart &gridPart_;
  const ProblemType &problem_;

};
#endif // #ifndef Stokes_MODEL_HH
