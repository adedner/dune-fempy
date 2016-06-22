#include <config.h>

// dune-fempy
#include <dune/fempy/python.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/pybind11/pybind11.h>

// dune-fem
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>

// dune-chns
#include "uzawascheme.hh"
#include "testingquantities.hh"

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

///////////////////////////////////////////////////////////////////////
// the wrapper for the Stokes scheme which we would like to expose to python
template <class StokesScheme>
struct StokesSchemeWrapper
{
  typedef StokesScheme BaseScheme;
  typedef typename StokesScheme::DiscreteFunctionType DiscreteFunction;
  typedef typename StokesScheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
  typedef typename StokesScheme::PressureDiscreteFunctionType PressureDiscreteFunction;
  typedef typename StokesScheme::ProblemType ProblemType;
  typedef typename StokesScheme::GridPartType GridPartType;
  typedef typename GridPartType::GridType HGridType;
  typedef typename StokesScheme::FullFunctionSpaceType FunctionSpaceType;

  StokesSchemeWrapper( GridPartType &gridPart, int problemNumber, double timestep ) :
    gridPart_( gridPart ),
    timeProvider_( gridPart_.grid() ),
    timestep_( timestep )
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
  StokesSchemeWrapper( StokesSchemeWrapper& ) = delete;
  StokesSchemeWrapper& operator=( const StokesSchemeWrapper& ) = delete;

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
    timeProvider_.next( timestep_ );
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
  const double timestepfactor_ = 0.29;
  double timestep_;
  const double factor_ = 0.5;
  const double viscosityActual_ = viscosity_*factor_;
  const double timestepStokes_ = 1/timestepfactor_;
  const double timestepBurgers_ = 1/( 1 - 2*timestepfactor_ );
  GridPartType &gridPart_;
  Dune::Fem::GridTimeProvider< HGridType > timeProvider_;
  ProblemType* problemPtr_ = 0;
  std::shared_ptr<StokesScheme> stokesScheme_;
};

namespace PyDune
{
  template <class Scheme>
  struct PythonExt
  {
    typedef Scheme StokesSchemeType;
    typedef typename Scheme::BaseScheme StokesScheme;
    typedef typename StokesScheme::ProblemType ProblemType;
    typedef typename StokesScheme::GridPartType GridPartType;
    typedef typename StokesScheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
    typedef typename StokesScheme::PressureDiscreteFunctionType PressureDiscreteFunction;

    static void updatevelocity( const StokesSchemeType *self, VelocityDiscreteFunction &velocity )
    {
      self->duneType()->updatevelocity( velocity );
    }
    static void updatepressure( const StokesSchemeType *self, PressureDiscreteFunction &pressure )
    {
      self->duneType()->updatepressure( pressure );
    }
    static void initialize(const StokesSchemeType *self)
    {
      self->duneType()->initialize();
    }
    static void preparestep1(const StokesSchemeType *self)
    {
      self->duneType()->preparestep1();
    }
    static void preparestep3(const StokesSchemeType *self)
    {
      self->duneType()->preparestep3();
    }
  };
  template< class StokesSchemeType, class... Args >
  bool addToPython (pybind11::class_< StokesSchemeType, Args... > &cls )
  {
    cls.def( "updatevelocity", &PythonExt<StokesSchemeType>::updatevelocity );
    cls.def( "updatepressure", &PythonExt<StokesSchemeType>::updatepressure );
    cls.def( "initialize",     &PythonExt<StokesSchemeType>::initialize );
    cls.def( "preparestep1",   &PythonExt<StokesSchemeType>::preparestep1 );
    cls.def( "preparestep3",   &PythonExt<StokesSchemeType>::preparestep3 );
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
      typedef StokesSchemeWrapper<Scheme> StokesSchemeType;
      typedef typename Scheme::GridPartType GridPartType;
      typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
      typedef typename Scheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
      typedef typename Scheme::PressureDiscreteFunctionType PressureDiscreteFunction;
      auto velo = detail::registerGridFunction< VelocityDiscreteFunction >( module, "VelocityDiscreteFunction" );
      auto pres = detail::registerGridFunction< PressureDiscreteFunction >( module, "PressureDiscreteFunction" );
      // export the scheme wrapper
      pybind11::class_< StokesSchemeType, std::shared_ptr<StokesSchemeType> > cls( module, "StokesScheme");
      cls.def( "__init__", [] ( StokesSchemeType &instance, GridPartType &gridPart, int modelNumber, double timestep ) {
          new( &instance ) StokesSchemeType( gridPart, modelNumber, timestep );
        }, pybind11::keep_alive< 1, 2 >() );
      cls.def( "velocity", [] (StokesSchemeType &scheme) -> VelocityDiscreteFunction& { return scheme.velocity(); },
            pybind11::return_value_policy::reference_internal );
      cls.def( "pressure", [] (StokesSchemeType &scheme) -> PressureDiscreteFunction& { return scheme.pressure(); },
            pybind11::return_value_policy::reference_internal );
      cls.def( "solve", &StokesSchemeType::solve );
      cls.def( "next", &StokesSchemeType::next );
      cls.def( "time", &StokesSchemeType::time );
      // add the user extensions
      PyDune::addToPython(cls);
    }
  }
}
