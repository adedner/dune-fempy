#include <config.h>

// dune-fempy
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/pybind11/pybind11.h>

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
  typedef NSBaseScheme<BurgersScheme> BaseType;
  typedef typename BurgersScheme::VelocitySpaceType VelocitySpaceType;
  typedef typename BurgersScheme::PressureSpaceType PressureSpaceType;
  typedef typename BurgersScheme::DiscreteFunctionType VelocityDiscreteFunction;
  typedef typename BurgersScheme::PressureDiscreteFunctionType PressureDiscreteFunction;
  typedef std::tuple<VelocitySpaceType&, PressureSpaceType&>
          SolutionSpaceType;
  typedef std::tuple<VelocityDiscreteFunction&, PressureDiscreteFunction&>
          SolutionType;

  using BaseType::model;
  BurgersSchemeWrapper( const SolutionSpaceType &spaces, double viscosity, int problemNumber, double timestep )
  : BaseType( std::get<0>(spaces).gridPart(), viscosity, problemNumber, timestep ),
    burgersScheme_ (std::get<0>(spaces),std::get<1>(spaces),model(),
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
    template< class Scheme >
    void registerScheme ( pybind11::module module )
    {
      typedef BurgersSchemeWrapper<Scheme> BurgersSchemeType;
      typedef typename BurgersSchemeType::SolutionSpaceType SolutionSpaceType;
      typedef typename BurgersSchemeType::SolutionType SolutionType;
      // export PRPScheme
      pybind11::class_< NSBaseScheme<Scheme> > clsBase( module, "NSBaseBScheme");
      pybind11::class_< BurgersSchemeType > cls( module, "Scheme", pybind11::base<NSBaseScheme<Scheme>>() );
      cls.def( "__init__", [] ( BurgersSchemeType &instance, const SolutionSpaceType &spaces,
                         double viscosity,
                         int problemNumber,
                         const std::string &name,
                         double timeStep ) {
          new( &instance ) BurgersSchemeType( spaces, viscosity, problemNumber, timeStep );
        }, pybind11::keep_alive< 1, 2 >(),
           pybind11::arg("spaces"),
           pybind11::arg("viscosity"),
           pybind11::arg("problemNumber"),
           pybind11::arg("name"),
           pybind11::arg("timeStep")
          );
      cls.def( "_solve",
          []( BurgersSchemeType &instance, const SolutionType &solution, bool assemble)
          { instance._solve(solution,assemble); } );
      cls.def( "_prepare", &BurgersSchemeType::_prepare );
      cls.def( "next", &BurgersSchemeType::next );
      cls.def( "time", &BurgersSchemeType::time );
    }
  }
}
