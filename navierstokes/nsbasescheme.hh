#ifndef NS_BASE_SCHEME
#define NS_BASE_SCHEME
#include "model.hh"
template <class Scheme>
struct NSBaseScheme
{
  typedef typename Scheme::ProblemType ProblemType;
  typedef typename Scheme::GridPartType GridPartType;
  typedef typename GridPartType::GridType HGridType;
  typedef typename Scheme::FullFunctionSpaceType FunctionSpaceType;
  typedef MyDiffusionModel<FunctionSpaceType,GridPartType> ModelType;

  NSBaseScheme( GridPartType &gridPart, double viscosity, int problemNumber, double timestep )
  : gridPart_( gridPart ),
    viscosity_(viscosity),
    timeProvider_( gridPart_.grid() ),
    timestep_( timestep )
  {
    timeProvider_.init( timestep_ );
    switch (problemNumber)
    {
      case 0:  problemPtr_ = new ChannelFlow<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
      case 1:  problemPtr_ = new VorticityFlow<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
      case 2:  problemPtr_ = new MovingPlate<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
      case 3:  problemPtr_ = new CouetteFlow<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
      case 4:  problemPtr_ = new KarmanVortexStreet<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
      default: problemPtr_ = new ChannelFlow<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
    }
    model_ = new ModelType(*problemPtr_); // , gridPart_);
  }
  ~NSBaseScheme() { delete problemPtr_; }
  const ModelType &model() const
  {
    return *model_;
  }
  void next()
  {
    timeProvider_.next( timestep_ );
  }
  double time()
  {
    return timeProvider_.time();
  }
  GridPartType &gridPart_;
  const double viscosity_; //  = 0.02;
  const double timestepfactor_ = 0.29289321881;
  const double factor_ = 0.58578643762;
  const double viscosityActual_ = viscosity_*factor_;
  const double timestepStokes_ = 1./timestepfactor_;
  const double timestepBurgers_ = 1./( 1. - 2.*timestepfactor_ );
  Dune::Fem::GridTimeProvider< HGridType > timeProvider_;
  double timestep_;
  ProblemType* problemPtr_ = 0;
  ModelType* model_;
};
#endif // NS_BASE_SCHEME
