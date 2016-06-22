#include <config.h>

// dune-fempy
#include <dune/fempy/python.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/pybind11/pybind11.h>

// dune-fem
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>

// dune-chns
#include "../../dune-chns/src/navier_stokes/prpscheme.hh"
#include "../../dune-chns/src/navier_stokes/testingquantities.hh"

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

//typedef PRPScheme< Dune::Fem::AdaptiveLeafGridPart< Dune::YaspGrid< 2, Dune::EquidistantCoordinates< double, 2 > > > > BurgersScheme;
typedef PRPScheme< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > > > BurgersScheme;

///////////////////////////////////////////////////////////////////////
// the wrapper for the Burgers scheme which we would like to expose to python
struct BurgersSchemeWrapper
{
  typedef typename BurgersScheme::DiscreteFunctionType DiscreteFunctionType;
  typedef Dune::GridSelector::GridType HGridType;
  typedef typename BurgersScheme::ProblemType ProblemType;
  typedef typename BurgersScheme::GridPartType GridPartType;
  typedef BurgersScheme::FullFunctionSpaceType FunctionSpaceType;

  BurgersSchemeWrapper( GridPartType &gridPart, int problemNumber ) :
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
      burgersScheme_ = std::make_shared<BurgersScheme>( gridPart_, *problemPtr_, timestepBurgers_, viscosityActual_ );
      timeProvider_.init( timestep_ );
    }
  ~BurgersSchemeWrapper() {std::cout << "BurgersSchemeWrapper destructor\n";
    delete problemPtr_;
  }
  BurgersSchemeWrapper(BurgersSchemeWrapper&) = delete;
  BurgersSchemeWrapper& operator=(const BurgersSchemeWrapper&) = delete;

  DiscreteFunctionType solution()
  {
    return duneType()->solution();
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
  std::shared_ptr<BurgersScheme> duneType() const
  {
    return burgersScheme_;
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
  std::shared_ptr<BurgersScheme> burgersScheme_;
};

namespace PyDune
{
  struct PythonExt2
  {
    typedef BurgersScheme::ProblemType ProblemType;
    typedef BurgersScheme::GridPartType GridPartType;

    // for velocity and pressure
    typedef BurgersScheme::VelocityDiscreteFunctionType VeloDF;
    typedef BurgersScheme::PressureDiscreteFunctionType PresDF;
    static PresDF pressure( const BurgersSchemeWrapper *self )
    {
      return self->duneType()->pressure();
    }
    static void updatevelocity( const BurgersSchemeWrapper *self, VeloDF velocity )
    {
      self->duneType()->updatevelocity( velocity );
    }
    static void updatepressure( const BurgersSchemeWrapper *self, PresDF pressure )
    {
      self->duneType()->updatepressure( pressure );
    }
    static void prepare(const BurgersSchemeWrapper *self)
    {
      self->duneType()->prepare();
    }
  };
  template< class... Args >
  bool addToPython (pybind11::class_< BurgersSchemeWrapper, Args... > &cls )
  {
    cls.def("updatevelocity", &PythonExt2::updatevelocity);
    cls.def("pressure", &PythonExt2::pressure, pybind11::return_value_policy::reference_internal);
    cls.def("updatepressure", &PythonExt2::updatepressure);
    cls.def("prepare", &PythonExt2::prepare);
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
      typedef typename BurgersScheme::GridPartType GridPartType;
      // export PRPScheme
      pybind11::class_< BurgersSchemeWrapper, std::shared_ptr<BurgersSchemeWrapper> > cls2( module, "BurgersScheme");
      cls2.def( "__init__", [] ( BurgersSchemeWrapper &instance, GridPartType &gridPart, int modelNumber ) {
          new( &instance ) BurgersSchemeWrapper( gridPart, modelNumber );
        }, pybind11::keep_alive< 1, 2 >() );
      cls2.def( "solve", &BurgersSchemeWrapper::solve );
      cls2.def( "solution", &BurgersSchemeWrapper::solution, pybind11::return_value_policy::reference_internal );
      cls2.def( "next", &BurgersSchemeWrapper::next );
      cls2.def( "time", &BurgersSchemeWrapper::time );
      // add the user extensions
      PyDune::addToPython(cls2);
    }
  }
}
