#ifndef DUNE_FEMPY_PYBIND11_GRIDFUNCTION_HH
#define DUNE_FEMPY_PYBIND11_GRIDFUNCTION_HH

#include <type_traits>

#include <dune/common/visibility.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    class HasLocalFunction;

  } // namespace Fem


  namespace FemPy
  {

    // getGridFunctionWrapper
    // ----------------------

    DUNE_EXPORT inline pybind11::object getGridFunctionWrapper ()
    {
      pybind11::object o = pybind11::module::import( "dune.ufl" ).attr( "GridFunction" );
      return o;
    }

  } // namespace FemPy

} // namespace Dune


#if 0
namespace pybind11
{

  namespace detail
  {

    // type_caster for dune-fem grid functions
    // ---------------------------------------

    template< class T >
    struct type_caster< T, std::enable_if_t< std::is_base_of< Dune::Fem::HasLocalFunction, T >::value > >
      : public type_caster_base< T >
    {
      typedef type_caster_base< T > Base;

      bool load ( handle src, bool convert )
      {
        bool ret;
        std::cout << "HANDLE A" << std::endl;
        if( isinstance( src, Dune::FemPy::getGridFunctionWrapper() ) )
        {
          std::cout << "HANDLE B-then" << std::endl;
          ret = Base::load( getattr( src, "__impl__" ), convert );
        }
        else
        {
          std::cout << "HANDLE B-else" << std::endl;
          ret = Base::load( src, convert );
        }
        std::cout << "HANDLE C" << std::endl;
        return ret;
      }

      template< class V >
      static handle cast ( V &&v, return_value_policy policy, handle parent )
      {
        std::cout << "CAST A" << std::endl;
        tuple args( 1 );
        std::cout << "CAST B" << std::endl;
        args[ 0 ] = Base::cast( std::forward< V >( v ), policy, parent );
        std::cout << "CAST C" << std::endl;
        handle ret = PyObject_Call( Dune::FemPy::getGridFunctionWrapper().ptr(), args.ptr(), nullptr );
        std::cout << "CAST D" << std::endl;
        return ret;
      }
    };

  } // namespace detail
} // namespace pybind11
#endif

#endif // #ifndef DUNE_FEMPY_PYBIND11_GRIDFUNCTION_HH
