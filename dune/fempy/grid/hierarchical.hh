#ifndef DUNE_FEMPY_GRID_HIERARCHICAL_HH
#define DUNE_FEMPY_GRID_HIERARCHICAL_HH

#include <cstddef>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/grid/common/rangegenerators.hh>

namespace Dune
{

  namespace FemPy
  {

    /*!
        @file
        @brief Contains C++ template classes for grids.
        \ingroup Grids
    */

    // readDGF
    // -------

    template< class Grid >
    static std::shared_ptr< Grid > readDGF ( const std::string &dgf )
    {
      GridPtr< Grid > gridPtr( dgf );
      gridPtr->loadBalance();
      return std::shared_ptr< Grid >( gridPtr.release() );
    }


    // HierarchicalGrid
    // ----------------

    template< class G >
    struct HierarchicalGrid
    {
      typedef G Grid;

      enum class Marker { Coarsen = -1, Keep = 0, Refine = 1 };

      typedef typename Grid::template Codim< 0 >::Entity Element;

      explicit HierarchicalGrid ( const std::string &dgf )
        : grid_( readDGF< Grid >( dgf ) )
      {}

      explicit HierarchicalGrid ( Grid *grid )
        : grid_( grid )
      {}

      template< class Marking >
      void mark ( Marking marking )
      {
        for( const Element &element : elements( grid_->leafGridView() ) )
        {
          Marker marker = marking( element );
          grid_->mark( static_cast< int >( marker ), element );
        }
      }

      void adapt ( )
      {
      }

      void globalRefine ( int level )
      {
        grid_->globalRefine( level );
      }

      void loadBalance ( )
      {
        grid_->loadBalance();
      }

      std::shared_ptr< Grid > grid () const { return grid_; }

    private:
      std::shared_ptr< Grid > grid_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_HIERARCHICAL_HH
