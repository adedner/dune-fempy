#include <config.h>

// dune-fempy
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/pybind11/pybind11.h>

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
  typedef NSBaseScheme<StokesScheme> BaseType;
  typedef typename StokesScheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
  typedef typename StokesScheme::PressureDiscreteFunctionType PressureDiscreteFunction;
  typedef typename StokesScheme::ProblemType ProblemType;
  typedef typename StokesScheme::GridPartType GridPartType;
  typedef std::tuple<VelocityDiscreteFunction&, PressureDiscreteFunction&>
          SolutionType;

  StokesSchemeWrapper( GridPartType &gridPart, int problemNumber, double timestep )
  : BaseType(gridPart,problemNumber,timestep)
  , stokesScheme_( BaseType::gridPart_, *BaseType::problemPtr_, BaseType::viscosityActual_, BaseType::timestepStokes_ )
  , solution_( stokesScheme_.velocity(), stokesScheme_.pressure() )
  {}
  ~StokesSchemeWrapper() {std::cout << "StokesSchemeWrapper destructor\n";
  }
  StokesSchemeWrapper( StokesSchemeWrapper& ) = delete;
  StokesSchemeWrapper& operator=( const StokesSchemeWrapper& ) = delete;

  SolutionType &solution()
  {
    return solution_;
  }
  void solve( bool assemble )
  {
    duneType().solve( assemble );
  }
  void update( const SolutionType &solution )
  {
    duneType().updatevelocity( std::get<0>(solution) );
    duneType().updatepressure( std::get<1>(solution) );
  }
  void initialize()
  {
    duneType().initialize();
  }
  void preparestep1()
  {
    duneType().preparestep1();
  }
  void preparestep3()
  {
    duneType().preparestep3();
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
  SolutionType solution_;
};

namespace Dune
{
  namespace FemPy
  {
    template< class Scheme >
    void registerScheme ( pybind11::module module )
    {
      typedef StokesSchemeWrapper<Scheme> StokesSchemeType;
      typedef typename Scheme::GridPartType GridPartType;
      typedef typename Scheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
      typedef typename Scheme::PressureDiscreteFunctionType PressureDiscreteFunction;
      auto velo = detail::registerGridFunction< VelocityDiscreteFunction >( module, "VelocityDiscreteFunction" );
      auto pres = detail::registerGridFunction< PressureDiscreteFunction >( module, "PressureDiscreteFunction" );
      // export the scheme wrapper
      pybind11::class_< NSBaseScheme<Scheme> > clsBase( module, "NSBaseSScheme");
      pybind11::class_< StokesSchemeType > cls( module, "StokesScheme", pybind11::base<NSBaseScheme<Scheme>>() );
      cls.def( "__init__", [] ( StokesSchemeType &instance, GridPartType &gridPart, int modelNumber, double timestep ) {
          new( &instance ) StokesSchemeType( gridPart, modelNumber, timestep );
        }, pybind11::keep_alive< 1, 2 >() );
      cls.def( "solution", &StokesSchemeType::solution,
            pybind11::return_value_policy::reference_internal );
      cls.def( "solve", &StokesSchemeType::solve );
      cls.def( "update", &StokesSchemeType::update );
      cls.def( "initialize", &StokesSchemeType::initialize );
      cls.def( "preparestep1", &StokesSchemeType::preparestep1 );
      cls.def( "preparestep3", &StokesSchemeType::preparestep3 );
      cls.def( "next", &StokesSchemeType::next );
      cls.def( "time", &StokesSchemeType::time );
    }
  }
}
