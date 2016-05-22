#ifndef DUNE_FEMPY_PY_GRID_RANGE_HH
#define DUNE_FEMPY_PY_GRID_RANGE_HH

#include <string>
#include <utility>

#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // PyGridViewRange
    // ---------------

    template< class GridView, int codim >
    struct PyGridViewRange
    {
      typedef typename GridView::template Codim< codim >::Iterator Iterator;

      PyGridViewRange ( const GridView &gridView, pybind11::object ref )
        : gridView_( gridView ), ref_( std::move( ref ) )
      {}

      Iterator begin () const { return gridView_.template begin< codim >(); }
      Iterator end () const { return gridView_.template end< codim >(); }

    private:
      const GridView &gridView_;
      pybind11::object ref_;
    };



    // PyGridViewIterator
    // ------------------

    template< class GridView, int codim >
    struct PyGridViewIterator
    {
      typedef PyGridViewRange< GridView, codim > Range;
      typedef typename GridView::template Codim< codim >::Entity Entity;

      PyGridViewIterator ( const Range &range ) : range_( range ), it_( range_.begin() ) {}

      Entity next ()
      {
        if( it_ == range_.end() )
          throw pybind11::stop_iteration();

        Entity entity = *it_;
        ++it_;
        return entity;
      }

    private:
      Range range_;
      typename Range::Iterator it_;
    };



    // registerPyGridViewRange
    // -----------------------

    template< class GridView, int codim >
    void registerPyGridViewRange ( pybind11::handle scope, const char *rgName )
    {
      typedef PyGridViewRange< GridView, codim > Range;
      typedef PyGridViewIterator< GridView, codim > Iterator;

      static const std::string itName = std::string( rgName ) + "Iterator";
      pybind11::class_< Iterator > itCls( scope, itName.c_str() );
      itCls.def( "__iter__", [] ( Iterator &it ) -> Iterator & { return it; } );
      itCls.def( "__next__", &Iterator::next );

      pybind11::class_< Range > rgCls( scope, rgName );
      rgCls.def( "__iter__", [] ( const Range &range ) { return Iterator( range ); } );
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_RANGE_HH
