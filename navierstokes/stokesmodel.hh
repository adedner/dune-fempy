#ifndef Stokes_MODEL_HH
#define Stokes_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

//#include "temporalprobleminterface.hh"
#include "navierstokes.hh"
#include "model.hh"


template <class Model>
struct StokesModel
{
  typedef Model ModelType;
  typedef typename ModelType::FunctionSpaceType FunctionSpaceType;
  typedef typename ModelType::GridPartType GridPartType;

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
  StokesModel( const ModelType &model,
               const GridPartType &gridPart,
               double mu, double nu)
    : model_(model),
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
    for (unsigned int i=0;i<GridPartType::dimensionworld;++i)
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
    bool isDirichlet = model_.isDirichletIntersection(inter, dirichletComponent);
    dirichletComponent[dimRange - 1] = false;

    return isDirichlet;
  }
protected:
  const ModelType &model_;
  double stab_;
  double mu_,nu_;

};

//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------

template< class Model >
struct StokesMainModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                                    Model::GridPartType::dimensionworld, Model::GridPartType::dimensionworld > ,
                                                    typename Model::GridPartType>
{
  typedef Model ModelType;
  typedef typename ModelType::GridPartType GridPartType;
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
  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  //! constructor
  StokesMainModel( const ModelType &model,
                   const GridPartType &gridPart,
                   double dt,
                   double mu, double nu,const bool implicit)
    : model_(model),
      dt_(dt),
      gridPart_(gridPart),
      zeroVelocityFunction_( "zero velocity", localZeroVelocity_, gridPart_ ),
      mu_(mu), nu_(nu),
      implicit_(implicit),
      dbc_(model_.dirichletBoundary()),
      nbc_(model_.neumanBoundary()),
      rhs_(model_.rightHandSide())
  {
    if (implicit==true)
    {
      timeStepFactor_=1.0;
      timeStepTheta_ =1.0;
    }
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
    for (unsigned int i=0;i<GridPartType::dimensionworld;++i)
      source[i] = (1./dt_)*nu_*value[i];
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
    return model_.hasDirichletBoundary() ;
  }
  bool hasNeumanBoundary () const
  {
    return model_.hasNeumanBoundary() ;
  }

  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<bool,dimRange> &dirichletComponent ) const
  {
    Dune::FieldVector<bool,dimRange+1> d;
    bool r = model_.isDirichletIntersection( inter, d );
    for (int i=0;i<dimRange;++i) dirichletComponent[i] = d[i];
    return r;
  }
  // return Fem :: Function for Dirichlet boundary values


  typedef GridPartType GP;
  struct LocalZeroVelocityFunction
  {
    typedef GP GridPartType;
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

public:
  const ZeroVelocityFunctionType &zeroVelocity() const
  {
    return zeroVelocityFunction_;
  }
  Dune::Fem::LocalFunctionAdapter< LocalVelocityExtractor<typename ModelType::DirichletBoundaryType> > dirichletBoundary(  ) const
  {
    return Dune::Fem::LocalFunctionAdapter< decltype(dbc_) >("right hand side", dbc_, gridPart_ );
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
  Dune::Fem::LocalFunctionAdapter< LocalVelocityExtractor<typename ModelType::RightHandSideType> > rightHandSideNext(  ) const
  {
    return Dune::Fem::LocalFunctionAdapter< decltype(rhs_) >("right hand side", rhs_, gridPart_ );
  }
  //------------------------------------------------------------------------------------------

private:
  double dt_;
  const ModelType &model_;
  const GridPartType &gridPart_;
  LocalZeroVelocityFunction localZeroVelocity_;
  ZeroVelocityFunctionType zeroVelocityFunction_;
  double mu_,nu_;
  double timeStepFactor_,timeStepTheta_;
  const bool implicit_;
  mutable LocalVelocityExtractor<typename ModelType::DirichletBoundaryType> dbc_;
  mutable LocalVelocityExtractor<typename ModelType::NeumanBoundaryType> nbc_;
  mutable LocalVelocityExtractor<typename ModelType::RightHandSideType> rhs_;
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
template< class Model >
struct StokesPrecondModel :
            public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,Model::GridPartType::dimensionworld,1> ,typename Model::GridPartType>
{
  typedef Model ModelType;
  typedef typename Model::GridPartType GridPartType;
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

public:
  //! constructor
  StokesPrecondModel( const ModelType& model, const GridPartType &gridPart,
                      double mu, double nu)
    : model_(model),
      gridPart_(gridPart),
      zeroVelocityFunction_( "zero velocity", localZeroVelocity_, gridPart_ ),
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
    Dune::FieldVector<bool,GridPartType::dimensionworld+1> d;
    return model_.isDirichletIntersection( inter, d );
    // for (int i=0;i<dimRange;++i) dirichletComponent[i] = d[i];
    // return r;
  }
protected:
  typedef GridPartType GP;
  struct LocalZeroVelocityFunction
  {
    typedef GP GridPartType;
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
  const ModelType &model_;
  const GridPartType &gridPart_;
  LocalZeroVelocityFunction localZeroVelocity_;
  ZeroVelocityFunctionType zeroVelocityFunction_;
  double mu_,nu_;
};

//---------------------------------------------------------------//
//---------------------------------------------------------------//
//---------------------------------------------------------------//


template< class Model >
struct StokesTransportModel : public DiffusionModelInterface<Dune::Fem::FunctionSpace<double,double,
                                            Model::GridPartType::dimensionworld, Model::GridPartType::dimensionworld >, typename Model::GridPartType>
{
  typedef Model ModelType;
  typedef typename Model::GridPartType GridPartType;
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
  StokesTransportModel(const ModelType &model, const GridPartType &gridPart)
  : gridPart_(gridPart), model_(model)
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
    Dune::FieldVector<bool,dimRange+1> d;
    bool r = model_.isDirichletIntersection( inter, d );
    for (int i=0;i<dimRange;++i) dirichletComponent[i] = d[i];
    return r;
  }
protected:

private:
  const GridPartType &gridPart_;
  const ModelType &model_;
};
#endif // #ifndef Stokes_MODEL_HH
