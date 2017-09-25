#ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
#define DUNE_FEMPY_PY_GRID_GRIDPART_HH

#include <cassert>

#include <map>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>

#include <dune/corepy/grid/hierarchical.hh>
#include <dune/corepy/grid/vtk.hh>

#include <dune/fempy/grid/gridpartadapter.hh>
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // GridModificationListener
      // ------------------------

      template< class Grid >
      class GridModificationListener final
        : public CorePy::GridModificationListener< Grid >
      {
        typedef Fem::DofManager< Grid > DofManager;

      public:
        GridModificationListener ( const Grid &grid )
          : dofManager_( DofManager::instance( grid ) )
        {}

        virtual void postModification ( const Grid &grid )
        {
          dofManager_.resize();
          dofManager_.compress();
        }

      private:
        DofManager &dofManager_;
      };



      // addGridModificationListener
      // ---------------------------

      template< class Grid >
      inline static void addGridModificationListener ( const Grid &grid )
      {
        typedef GridModificationListener< Grid > Listener;
        for( const auto &listener : CorePy::detail::gridModificationListeners( grid ) )
        {
          if( dynamic_cast< Listener * >( listener.second ) )
            return;
        }

        pybind11::handle nurse = pybind11::detail::get_object_handle( &grid, pybind11::detail::get_type_info( typeid( Grid ) ) );
        CorePy::detail::addGridModificationListener( grid, new Listener( grid ), nurse );
      }



      // GridPartConverter
      // -----------------

      template< class GV >
      struct GridPartConverter
      {
        typedef GV GridView;
        typedef GridPartAdapter< GV > GridPart;

        GridPart &operator() ( pybind11::handle gridView )
        {
          auto result = instances_.emplace( gridView.ptr(), nullptr );
          auto pos = result.first;
          if( result.second )
          {
            GridView view = gridView.template cast< GridView >();

            // add grid modification listener (if not registered)
            addGridModificationListener( view.grid() );

            // create new gridpart object
            pos->second = new GridPart( view );

            // create Python guard object, removing the grid part once the grid view dies
            pybind11::cpp_function remove_gridpart( [ this, pos ] ( pybind11::handle weakref ) {
                delete pos->second;
                instances_.erase( pos );
                weakref.dec_ref();
              } );
            pybind11::weakref weakref( gridView, remove_gridpart );
            weakref.release();
          }
          assert( pos->second );
          return *pos->second;
        }

      private:
        std::map< PyObject *, GridPart * > instances_;
      };


      template< class GP >
      struct GridPartConverter< Dune::GridView< Fem::GridPart2GridViewTraits< GP > > >
      {
        typedef GP GridPart;
        typedef Dune::GridView< Fem::GridPart2GridViewTraits< GP > > GridView;

        GridPart &operator() ( pybind11::handle gridView )
        {
          return const_cast< GridPart & >( gridView.template cast< GridView >().impl().gridPart() );
        }
      };


      template< class GP >
      struct GridPartConverter< Fem::GridPart2GridViewImpl< GP > >
      {
        typedef GP GridPart;
        typedef Fem::GridPart2GridViewImpl< GP > GridView;

        GridPart &operator() ( pybind11::handle gridView )
        {
          return const_cast< GridPart & >( gridView.template cast< GridView >().gridPart() );
        }
      };



      // gridPartConverter
      // -----------------

      template< class GridView >
      inline GridPartConverter< GridView > &gridPartConverter ()
      {
        static GridPartConverter< GridView > converter;
        return converter;
      }

    } // namespace detail



    // GridPart
    // --------

    template< class GridView >
    using GridPart = typename detail::GridPartConverter< GridView >::GridPart;



    // gridPart
    // --------

    template< class GridView >
    inline static GridPart< GridView > &gridPart ( pybind11::handle gridView )
    {
      return detail::gridPartConverter< GridView >()( std::move( gridView ) );
    }




    // GridPartDeleter
    // ---------------

    template< class GridView >
    struct GridPartDeleter;

    template< class GridPart >
    struct GridPartDeleter< Dune::GridView< Fem::GridPart2GridViewTraits< GridPart > > >
    {
      void operator() ( Dune::GridView< Fem::GridPart2GridViewTraits< GridPart > > *gridView )
      {
        delete &gridView.template cast< GridView >().impl().gridPart();
        delete gridView;
      }
    };

    template< class GridPart >
    struct GridPartDeleter< Fem::GridPart2GridViewImpl< GridPart > >
    {
      void operator() ( Fem::GridPart2GridViewImpl< GridPart > *gridView )
      {
        delete &gridView.template cast< GridView >().gridPart();
        delete gridView;
      }
    };



    // GridPartPtr
    // -----------

    template< class GridView >
    using GridPartPtr = std::unique_ptr< GridView, GridPartDeleter< GridView > >;



    // makeGridPart
    // ------------

    template< class GridView, class... Args >
    inline static GridPartPtr< GridView > makeGridPart ( Args &&... args )
    {
      GridPart< GridView > *gridPart = new GridPart< GridView >( std::forward< Args >( args )... );
      return GridPartPtr< GridView >( new GridView( static_cast< GridView >( *gridPart ) ) );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
