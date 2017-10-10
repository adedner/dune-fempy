#ifndef DUNE_FEMPY_FUNCTION_UTILITY_HH
#define DUNE_FEMPY_FUNCTION_UTILITY_HH

#include <dune/common/typeutilities.hh>

namespace Dune
{

  namespace FemPy
  {

    // order
    // -----

    template< class GridFunction >
    inline static auto order ( const GridFunction &gridFunction, PriorityTag< 2 > )
      -> decltype( gridFunction.space().order() )
    {
      return gridFunction.space().order();
    }

    template< class GridFunction >
    inline static auto order ( const GridFunction &gridFunction, PriorityTag< 1 > )
      -> decltype( gridFunction.order() )
    {
      return gridFunction.order();
    }

    template< class GridFunction >
    inline static int order ( const GridFunction &gridFunction, PriorityTag< 0 > )
    {
      return std::numeric_limits< int >::max();
    }

    template< class GridFunction >
    inline static int order ( const GridFunction &gridFunction )
    {
      return order( gridFunction, PriorityTag< 42 >() );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_UTILITY_HH
