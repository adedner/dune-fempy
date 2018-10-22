#include <config.h>

// dune-python
#include <dune/python/pybind11/pybind11.h>

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
  typedef std::tuple<VelocitySpaceType&, PressureSpaceType&>
          DiscreteFunctionSpaceType;
  typedef std::tuple<VelocityDiscreteFunction&, PressureDiscreteFunction&>
          SolutionType;

  StokesSchemeWrapper( const DiscreteFunctionSpaceType &spaces, const ModelType &model, double viscosity, double timestep )
  : BaseType( viscosity, timestep )
  , stokesScheme_( std::get<0>(spaces), std::get<1>(spaces), model,
      BaseType::timestep_, BaseType::viscosityActual_, BaseType::timestepStokes_ )
  {}
  ~StokesSchemeWrapper() {std::cout << "StokesSchemeWrapper destructor\n";
  }
  StokesSchemeWrapper( StokesSchemeWrapper& ) = delete;
  StokesSchemeWrapper& operator=( const StokesSchemeWrapper& ) = delete;

  void _solve( const SolutionType &target, bool assemble )
  {
    duneType().solve( std::get<0>(target), std::get<1>(target), assemble );
  }
  void initialize( const SolutionType &solution )
  {
    duneType().initialize( std::get<0>(solution), std::get<1>(solution) );
  }
  void _prepare( const SolutionType &solution )
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

namespace Dune
{
  namespace FemPy
  {
    namespace detail
    {
      template< class Scheme, class Cls >
      void registerScheme ( pybind11::module module, Cls &cls, std::false_type )
      {
        typedef Scheme StokesSchemeType;
        typedef typename StokesSchemeType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        typedef typename StokesSchemeType::SolutionType SolutionType;
        typedef typename StokesSchemeType::ModelType ModelType;
        // export the scheme wrapper
        //pybind11::class_< NSBaseScheme<Scheme> > clsBase( module, "NSBaseSScheme");
        //pybind11::class_< StokesSchemeType > cls2( module, "Scheme", pybind11::base<NSBaseScheme<Scheme>>() );
        cls.def( "__init__", [] ( StokesSchemeType &instance, const DiscreteFunctionSpaceType &spaces,
                           const ModelType &model,
                           const std::string &name,
                           double viscosity,
                           double timeStep ) {
            new( &instance ) StokesSchemeType( spaces, model, viscosity, timeStep );
          }, pybind11::keep_alive< 1, 2 >(),
             pybind11::arg("spaces"),
             pybind11::arg("model"),
             pybind11::arg("name"),
             pybind11::arg("viscosity"),
             pybind11::arg("timeStep")
            );
        cls.def( "_solve",
            []( StokesSchemeType &instance, const SolutionType &solution, bool assemble)
            { instance._solve(solution,assemble); } );
        cls.def( "initialize", &StokesSchemeType::initialize );
        cls.def( "_prepare", &StokesSchemeType::_prepare );
      }
    }
  }
}