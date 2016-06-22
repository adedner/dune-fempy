#include <config.h>

// dune-fempy
#include <dune/fempy/python.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/pybind11/pybind11.h>

// dune-fem
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>

// dune-chns
#include "../../dune-chns/src/navier_stokes/uzawascheme.hh"
#include "../../dune-chns/src/navier_stokes/testingquantities.hh"

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

//typedef UzawaScheme< Dune::Fem::AdaptiveLeafGridPart< Dune::YaspGrid< 2, Dune::EquidistantCoordinates< double, 2 > > > > StokesScheme;
typedef UzawaScheme< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > > > StokesScheme;

///////////////////////////////////////////////////////////////////////
// the wrapper for the Stokes scheme which we would like to expose to python
struct StokesSchemeWrapper
{
  typedef typename StokesScheme::DiscreteFunctionType DiscreteFunction;
  typedef typename StokesScheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
  typedef typename StokesScheme::PressureDiscreteFunctionType PressureDiscreteFunction;
  typedef Dune::GridSelector::GridType HGridType;
  typedef typename StokesScheme::ProblemType ProblemType;
  typedef typename StokesScheme::GridPartType GridPartType;
  typedef StokesScheme::FullFunctionSpaceType FunctionSpaceType;

  StokesSchemeWrapper( GridPartType &gridPart, int problemNumber ) :
    gridPart_( gridPart ),
    timeProvider_( gridPart_.grid() )
    {
      switch (problemNumber)
      {
        case 0: problemPtr_ = new ChannelFlow<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
        case 1: problemPtr_ = new VorticityFlow<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
        case 2: problemPtr_ = new MovingPlate<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
        case 3: problemPtr_ = new CouetteFlow<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
        case 4: problemPtr_ = new KarmanVortexStreet<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
        default: problemPtr_ = new ChannelFlow<FunctionSpaceType>( timeProvider_, viscosity_, timestepStokes_, timestepBurgers_ ); break;
      }
      stokesScheme_ = std::make_shared<StokesScheme>( gridPart_, *problemPtr_, viscosityActual_, timestepBurgers_ );
      timeProvider_.init( timestep_ );
    }
  ~StokesSchemeWrapper() {std::cout << "StokesSchemeWrapper destructor\n";
    delete problemPtr_;
  }
  StokesSchemeWrapper(StokesSchemeWrapper&) = delete;
  StokesSchemeWrapper& operator=(const StokesSchemeWrapper&) = delete;

  DiscreteFunction &solution()
  {
    return duneType()->solution();
  }
  VelocityDiscreteFunction &velocity()
  {
    return duneType()->velocity();
  }
  PressureDiscreteFunction &pressure()
  {
    return duneType()->pressure();
  }
  void solve( bool assemble )
  {
    duneType()->solve( assemble );
  }
  void next()
  {
    timeProvider_.next(timestep_);
  }
  void time()
  {
    std::cout << timeProvider_.time() << std::endl;
  }
  std::shared_ptr<StokesScheme> duneType() const
  {
    return stokesScheme_;
  }
  protected:
  const double viscosity_ = 0.01;
  const double timestep_ = 0.29;
  const double factor_ = 0.5;
  const double viscosityActual_ = viscosity_*factor_;
  const double timestepStokes_ = 1/timestep_;
  const double timestepBurgers_ = 1/(1-2*timestep_);
  GridPartType &gridPart_;
  Dune::Fem::GridTimeProvider< HGridType > timeProvider_;
  ProblemType* problemPtr_ = 0;
  std::shared_ptr<StokesScheme> stokesScheme_;
};

namespace PyDune
{
  struct PythonExt
  {
    typedef StokesScheme::ProblemType ProblemType;
    typedef StokesScheme::GridPartType GridPartType;
    typedef StokesScheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
    typedef StokesScheme::PressureDiscreteFunctionType PressureDiscreteFunction;

    static void updatevelocity( const StokesSchemeWrapper *self, VelocityDiscreteFunction &velocity )
    {
      self->duneType()->updatevelocity( velocity );
    }
    static void updatepressure( const StokesSchemeWrapper *self, PressureDiscreteFunction &pressure )
    {
      self->duneType()->updatepressure( pressure );
    }
    static void initialize(const StokesSchemeWrapper *self)
    {
      self->duneType()->initialize();
    }
    static void preparestep1(const StokesSchemeWrapper *self)
    {
      self->duneType()->preparestep1();
    }
    static void preparestep3(const StokesSchemeWrapper *self)
    {
      self->duneType()->preparestep3();
    }
  };
  template< class... Args >
  bool addToPython (pybind11::class_< StokesSchemeWrapper, Args... > &cls )
  {
    cls.def( "updatevelocity", &PythonExt::updatevelocity );
    cls.def( "updatepressure", &PythonExt::updatepressure );
    cls.def( "initialize", &PythonExt::initialize );
    cls.def( "preparestep1", &PythonExt::preparestep1 );
    cls.def( "preparestep3", &PythonExt::preparestep3 );
    return false;
  }
}

namespace Dune
{

  namespace FemPy
  {
    template< class Scheme >
    void registerScheme ( pybind11::module module )
    {
      typedef typename StokesScheme::GridPartType GridPartType;
      typedef typename StokesScheme::DiscreteFunctionType DiscreteFunction;
      typedef typename StokesScheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
      typedef typename StokesScheme::PressureDiscreteFunctionType PressureDiscreteFunction;
      // auto sol = detail::registerGridFunction< DiscreteFunction >( module, "DiscreteFunction" );
      auto velo = detail::registerGridFunction< VelocityDiscreteFunction >( module, "VelocityDiscreteFunction" );
      auto pres = detail::registerGridFunction< PressureDiscreteFunction >( module, "PressureDiscreteFunction" );
      // export the scheme wrapper
      pybind11::class_< StokesSchemeWrapper, std::shared_ptr<StokesSchemeWrapper> > cls( module, "StokesScheme");
      cls.def( "__init__", [] ( StokesSchemeWrapper &instance, GridPartType &gridPart, int modelNumber ) {
          new( &instance ) StokesSchemeWrapper( gridPart, modelNumber );
        }, pybind11::keep_alive< 1, 2 >() );
      cls.def( "velocity", [] (StokesSchemeWrapper &scheme) -> VelocityDiscreteFunction& { return scheme.velocity(); },
            pybind11::return_value_policy::reference_internal );
      cls.def( "pressure", [] (StokesSchemeWrapper &scheme) -> PressureDiscreteFunction& { return scheme.pressure(); },
            pybind11::return_value_policy::reference_internal );
      cls.def( "solve", &StokesSchemeWrapper::solve );
      //cls.def( "solution", [] (StokesSchemeWrapper &scheme) -> DiscreteFunction& { return scheme.solution(); },
      //      pybind11::return_value_policy::reference_internal );
      cls.def( "next", &StokesSchemeWrapper::next );
      cls.def( "time", &StokesSchemeWrapper::time );
      // add the user extensions
      PyDune::addToPython(cls);
    }
  }
}
