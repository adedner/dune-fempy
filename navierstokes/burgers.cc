#include <config.h>

// dune-fempy
#include <dune/fempy/python.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/pybind11/pybind11.h>

// dune-fem
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>

// dune-chns
#include "prpscheme.hh"
#include "testingquantities.hh"

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

///////////////////////////////////////////////////////////////////////
// the wrapper for the Burgers scheme which we would like to expose to python
template <class BurgersScheme>
struct BurgersSchemeWrapper
{
  typedef BurgersScheme BaseScheme;
  typedef typename BurgersScheme::DiscreteFunctionType DiscreteFunction;
  typedef typename BurgersScheme::ProblemType ProblemType;
  typedef typename BurgersScheme::GridPartType GridPartType;
  typedef typename GridPartType::GridType HGridType;
  typedef typename BurgersScheme::FullFunctionSpaceType FunctionSpaceType;

  BurgersSchemeWrapper( GridPartType &gridPart, int problemNumber, double timestep ) :
    gridPart_( gridPart ),
    timeProvider_( gridPart_.grid() ),
    timestep_( timestep)
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

  DiscreteFunction &solution()
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
  double time()
  {
    return timeProvider_.time();
  }
  std::shared_ptr<BurgersScheme> duneType() const
  {
    return burgersScheme_;
  }
  protected:
  const double viscosity_ = 0.0001;
  const double timestepfactor_ = 0.29;
  double timestep_;
  const double factor_ = 0.585756;
  const double viscosityActual_ = viscosity_*factor_;
  const double timestepStokes_ = 1./timestepfactor_;
  const double timestepBurgers_ = 1./( 1. - 2.*timestepfactor_ );
  GridPartType &gridPart_;
  Dune::Fem::GridTimeProvider< HGridType > timeProvider_;
  ProblemType* problemPtr_ = 0;
  std::shared_ptr<BurgersScheme> burgersScheme_;
};

namespace PyDune
{
  template <class Scheme>
  struct PythonExt2
  {
    typedef Scheme BurgersSchemeType;
    typedef typename Scheme::BaseScheme BurgersScheme;
    typedef typename BurgersScheme::ProblemType ProblemType;
    typedef typename BurgersScheme::GridPartType GridPartType;
    typedef typename BurgersScheme::VelocityDiscreteFunctionType VelocityDiscreteFunction;
    typedef typename BurgersScheme::PressureDiscreteFunctionType PressureDiscreteFunction;

    static void updatevelocity( const BurgersSchemeType *self, VelocityDiscreteFunction velocity )
    {
      self->duneType()->updatevelocity( velocity );
    }
    static void updatepressure( const BurgersSchemeType *self, PressureDiscreteFunction pressure )
    {
      self->duneType()->updatepressure( pressure );
    }
    static void prepare(const BurgersSchemeType *self)
    {
      self->duneType()->prepare();
    }
  };
  template< class BurgersSchemeType, class... Args >
  bool addToPython (pybind11::class_< BurgersSchemeType, Args... > &cls )
  {
    cls.def("updatevelocity", &PythonExt2<BurgersSchemeType>::updatevelocity);
    cls.def("updatepressure", &PythonExt2<BurgersSchemeType>::updatepressure);
    cls.def("prepare", &PythonExt2<BurgersSchemeType>::prepare);
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
      typedef BurgersSchemeWrapper<Scheme> BurgersSchemeType;
      typedef typename Scheme::GridPartType GridPartType;
      typedef typename Scheme::DiscreteFunctionType SolutionFunction;
      // export PRPScheme
      pybind11::class_< BurgersSchemeType, std::shared_ptr<BurgersSchemeType> > cls2( module, "BurgersScheme");
      cls2.def( "__init__", [] ( BurgersSchemeType &instance, GridPartType &gridPart, int modelNumber, double timestep ) {
          new( &instance ) BurgersSchemeType( gridPart, modelNumber, timestep );
        }, pybind11::keep_alive< 1, 2 >() );
      cls2.def( "solve", &BurgersSchemeType::solve );
      cls2.def( "solution", [] (BurgersSchemeType &scheme) -> SolutionFunction& { return scheme.solution(); },
            pybind11::return_value_policy::reference_internal );
      cls2.def( "next", &BurgersSchemeType::next );
      cls2.def( "time", &BurgersSchemeType::time );
      // add the user extensions
      PyDune::addToPython(cls2);
    }
  }
}
