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

  NSBaseScheme(double viscosity, double timestep )
  : viscosity_(viscosity),
    timestep_( timestep )
  {
  }
  const double viscosity_; //  = 0.02;
  const double timestepfactor_ = 0.29289321881;
  const double factor_ = 0.58578643762;
  const double viscosityActual_ = viscosity_*factor_;
  const double timestepStokes_ = 1./timestepfactor_;
  const double timestepBurgers_ = 1./( 1. - 2.*timestepfactor_ );
  double timestep_;
};
#endif // NS_BASE_SCHEME
