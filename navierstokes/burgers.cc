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
#include "nsbasescheme.hh"
#include "testingquantities.hh"

///////////////////////////////////////////////////////////////////////
// the wrapper for the Burgers scheme which we would like to expose to python
template <class BurgersScheme>
struct BurgersSchemeWrapper : NSBaseScheme<BurgersScheme>
{
  typedef NSBaseScheme<BurgersScheme> BaseType;
  typedef BurgersScheme BaseScheme;
  typedef typename BurgersScheme::DiscreteFunctionType VelocityDiscreteFunction;
  typedef typename BurgersScheme::PressureDiscreteFunctionType PressureDiscreteFunction;
  typedef typename BurgersScheme::ProblemType ProblemType;
  typedef typename BurgersScheme::GridPartType GridPartType;
  typedef std::tuple<VelocityDiscreteFunction&, PressureDiscreteFunction&>
          SolutionType;

  BurgersSchemeWrapper( GridPartType &gridPart, int problemNumber, double timestep )
  : BaseType( gridPart, problemNumber, timestep ),
    burgersScheme_ (BaseType::gridPart_, *BaseType::problemPtr_, BaseType::timestepBurgers_, BaseType::viscosityActual_ )
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
  void solve( bool assemble )
  {
    duneType().solve( assemble );
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
    typedef typename Scheme::SolutionType SolutionType;

    static void update( BurgersSchemeType &self, const SolutionType &solution )
    {
      self.duneType().updatevelocity( std::get<0>(solution) );
    }
    static void prepare( BurgersSchemeType &self)
    {
      self.duneType().prepare();
    }
  };
  template< class BurgersSchemeType, class... Args >
  bool addToPython (pybind11::class_< BurgersSchemeType, Args... > &cls )
  {
    cls.def("update", &PythonExt2<BurgersSchemeType>::update);
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
      pybind11::class_< NSBaseScheme<Scheme> > clsBase( module, "NSBaseBScheme");
      pybind11::class_< BurgersSchemeType > cls2( module, "BurgersScheme",
          pybind11::base<NSBaseScheme<Scheme>>() );
      cls2.def( "__init__", [] ( BurgersSchemeType &instance, GridPartType &gridPart, int modelNumber, double timestep ) {
          new( &instance ) BurgersSchemeType( gridPart, modelNumber, timestep );
        }, pybind11::keep_alive< 1, 2 >() );
      cls2.def( "solve", &BurgersSchemeType::solve );
      cls2.def( "solution", &BurgersSchemeType::solution,
            pybind11::return_value_policy::reference_internal );
      cls2.def( "next", &BurgersSchemeType::next );
      cls2.def( "time", &BurgersSchemeType::time );
      // add the user extensions
      PyDune::addToPython(cls2);
    }
  }
}
