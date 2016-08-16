#include <config.h>

// dune-corepy
#include <dune/corepy/pybind11/pybind11.h>

// dune-fempy
#include <dune/fempy/py/grid/function.hh>

// dune-fem
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>

// dune-chns
#include "prpscheme.hh"
#include "nsbasescheme.hh"
#include "testingquantities.hh"

///////////////////////////////////////////////////////////////////////
// the wrapper for the Burgers scheme which we would like to expose to python
template <class BurgersScheme>
struct BurgersSchemeWrapper : NSBaseScheme<BurgersScheme>
{
  typedef typename BurgersScheme::AdditionalModelType ModelType;
  typedef NSBaseScheme<BurgersScheme> BaseType;
  typedef typename BurgersScheme::VelocitySpaceType VelocitySpaceType;
  typedef typename BurgersScheme::PressureSpaceType PressureSpaceType;
  typedef typename BurgersScheme::DiscreteFunctionType VelocityDiscreteFunction;
  typedef typename BurgersScheme::PressureDiscreteFunctionType PressureDiscreteFunction;
  typedef std::tuple<VelocitySpaceType&, PressureSpaceType&>
          DiscreteFunctionSpaceType;
  typedef std::tuple<VelocityDiscreteFunction&, PressureDiscreteFunction&>
          SolutionType;

  BurgersSchemeWrapper( const DiscreteFunctionSpaceType &spaces, const ModelType &model, double viscosity, double timestep )
  : BaseType( viscosity, timestep ),
    burgersScheme_ (std::get<0>(spaces),std::get<1>(spaces),model,
         BaseType::timestep_, BaseType::viscosityActual_, BaseType::timestepBurgers_ )
  {
  }
  ~BurgersSchemeWrapper() {std::cout << "BurgersSchemeWrapper destructor\n";
  }
  BurgersSchemeWrapper(BurgersSchemeWrapper&) = delete;
  BurgersSchemeWrapper& operator=(const BurgersSchemeWrapper&) = delete;

  void _solve( const SolutionType &target, bool assemble )
  {
    duneType().solve( std::get<0>(target), std::get<1>(target), assemble );
  }
  void _prepare( const SolutionType &solution )
  {
    duneType().prepare( std::get<0>(solution) );
  }
  const BurgersScheme& duneType() const
  {
    return burgersScheme_;
  }
  BurgersScheme& duneType()
  {
    return burgersScheme_;
  }
  protected:
  BurgersScheme burgersScheme_;
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
        typedef Scheme BurgersSchemeType;
        typedef typename BurgersSchemeType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        typedef typename BurgersSchemeType::SolutionType SolutionType;
        typedef typename BurgersSchemeType::ModelType ModelType;
        cls.def( "__init__", [] ( BurgersSchemeType &instance, const DiscreteFunctionSpaceType &spaces,
                           const ModelType &model,
                           const std::string &name,
                           double viscosity,
                           double timeStep ) {
            new( &instance ) BurgersSchemeType( spaces, model, viscosity, timeStep );
          }, pybind11::keep_alive< 1, 2 >(),
             pybind11::arg("spaces"),
             pybind11::arg("model"),
             pybind11::arg("name"),
             pybind11::arg("viscosity"),
             pybind11::arg("timeStep")
            );
        cls.def( "_solve",
            []( BurgersSchemeType &instance, const SolutionType &solution, bool assemble)
            { instance._solve(solution,assemble); } );
        cls.def( "_prepare", &BurgersSchemeType::_prepare );
      }
    }
  }
}
