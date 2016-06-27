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
  typedef BurgersScheme BaseScheme;
  typedef typename BurgersScheme::VelocityDiscreteFunctionSpaceType VelocityDiscreteFunctionSpaceType;
  typedef typename BurgersScheme::PressureDiscreteFunctionSpaceType PressureDiscreteFunctionSpaceType;
  typedef typename BurgersScheme::DiscreteFunctionType VelocityDiscreteFunction;
  typedef typename BurgersScheme::PressureDiscreteFunctionType PressureDiscreteFunction;
  typedef typename BurgersScheme::ProblemType ProblemType;
  typedef typename BurgersScheme::GridPartType GridPartType;
  typedef std::tuple<VelocityDiscreteFunctionSpaceType&, PressureDiscreteFunctionSpaceType&>
          SolutionSpaceType;
  typedef std::tuple<VelocityDiscreteFunction&, PressureDiscreteFunction&>
          SolutionType;

  BurgersSchemeWrapper( const SolutionSpaceType &spaces, int problemNumber, double timestep )
  : BaseType( std::get<0>(spaces).gridPart(), problemNumber, timestep ),
    burgersScheme_ (std::get<0>(spaces),std::get<1>(spaces), *BaseType::problemPtr_, BaseType::timestepBurgers_, BaseType::viscosityActual_ )
  , solution_( burgersScheme_.solution(), burgersScheme_.pressure() )
  {
  }
  ~BurgersSchemeWrapper() {std::cout << "BurgersSchemeWrapper destructor\n";
  }
  BurgersSchemeWrapper(BurgersSchemeWrapper&) = delete;
  BurgersSchemeWrapper& operator=(const BurgersSchemeWrapper&) = delete;

  SolutionType &solution()
  {
    return solution_;
  }
  void _solve( const SolutionType &target, bool assemble )
  {
    duneType().solve( assemble );
  }
  void _prepare( const SolutionType &solution )
  {
    duneType().updatevelocity( std::get<0>(solution) );
    duneType().prepare();
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
  SolutionType solution_;
};

namespace Dune
{
  namespace FemPy
  {
    template< class Scheme >
    void registerScheme ( pybind11::module module )
    {
      typedef BurgersSchemeWrapper<Scheme> BurgersSchemeType;
      typedef typename Scheme::GridPartType GridPartType;
      typedef typename BurgersSchemeType::SolutionSpaceType SolutionSpaceType;
      typedef typename Scheme::DiscreteFunctionType SolutionFunction;
      // export PRPScheme
      pybind11::class_< NSBaseScheme<Scheme> > clsBase( module, "NSBaseBScheme");
      pybind11::class_< BurgersSchemeType > cls( module, "Scheme", pybind11::base<NSBaseScheme<Scheme>>() );
      cls.def( "__init__", [] ( BurgersSchemeType &instance, const SolutionSpaceType &spaces,
                                int modelNumber, double timestep ) {
          new( &instance ) BurgersSchemeType( spaces, modelNumber, timestep );
        }, pybind11::keep_alive< 1, 2 >() );
      cls.def( "_solve", &BurgersSchemeType::_solve );
      cls.def( "solution", &BurgersSchemeType::solution,
            pybind11::return_value_policy::reference_internal );
      cls.def( "_prepare", &BurgersSchemeType::_prepare );
      cls.def( "next", &BurgersSchemeType::next );
      cls.def( "time", &BurgersSchemeType::time );
    }
  }
}
