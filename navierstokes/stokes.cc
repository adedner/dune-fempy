#include <config.h>

// dune-corepy
#include <dune/corepy/pybind11/pybind11.h>

// dune-fempy
#include <dune/fempy/py/grid/function.hh>

// dune-fem
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>

// dune-chns
#include "uzawascheme.hh"
#include "nsbasescheme.hh"
#include "testingquantities.hh"

///////////////////////////////////////////////////////////////////////
// the wrapper for the Stokes scheme which we would like to expose to python
template <class StokesScheme>
struct StokesSchemeWrapper : public NSBaseScheme<StokesScheme>
{
  typedef typename StokesScheme::AdditionalModelType ModelType;
  typedef NSBaseScheme<StokesScheme> BaseType;
  typedef typename StokesScheme::VelocitySpaceType VelocitySpaceType;
  typedef typename StokesScheme::PressureSpaceType PressureSpaceType;
  typedef typename StokesScheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
  typedef typename StokesScheme::PressureDiscreteFunctionType PressureDiscreteFunction;
  typedef std::tuple<VelocitySpaceType&, PressureSpaceType&> DiscreteFunctionSpaceType;
  typedef std::tuple<VelocityDiscreteFunction&, PressureDiscreteFunction&>
          DiscreteFunctionType;
  typedef typename StokesScheme::GridPartType GridPartType;

  StokesSchemeWrapper( const DiscreteFunctionSpaceType &spaces, const ModelType &model, double viscosity, double timestep )
  : BaseType( viscosity, timestep )
  , stokesScheme_( std::get<0>(spaces), std::get<1>(spaces), model,
      BaseType::timestep_, BaseType::viscosityActual_, BaseType::timestepStokes_ )
  {}
  ~StokesSchemeWrapper() {std::cout << "StokesSchemeWrapper destructor\n";
  }
  StokesSchemeWrapper( StokesSchemeWrapper& ) = delete;
  StokesSchemeWrapper& operator=( const StokesSchemeWrapper& ) = delete;

  void _solve( const DiscreteFunctionType &target, bool assemble )
  {
    duneType().solve( std::get<0>(target), std::get<1>(target), assemble );
  }
  void initialize( const DiscreteFunctionType &solution )
  {
    duneType().initialize( std::get<0>(solution), std::get<1>(solution) );
  }
  void _prepare( const DiscreteFunctionType &solution )
  {
    duneType().prepare( std::get<0>(solution) );
  }
  StokesScheme& duneType()
  {
    return stokesScheme_;
  }
  const StokesScheme& duneType() const
  {
    return stokesScheme_;
  }
  protected:
  StokesScheme stokesScheme_;
};
