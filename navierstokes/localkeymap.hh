#ifndef SPACE_LOCALKEYMAP_HH
#define SPACE_LOCALKEYMAP_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/space/mapper/code.hh>
#include <dune/fem/space/mapper/localkey.hh>

namespace Dune
{

  namespace Fem
  {

    template< int dim >
    struct BubbleElementLocalKeyMap
    {
      //! [Constructor of LocalKey tripple]
      BubbleElementLocalKeyMap ()
      {
        for( int i = 0; i <= dim; ++i )
          map_.emplace_back( i, dim, 0 );
        map_.emplace_back( 0, 0, 0 );
      }
      //! [Constructor of LocalKey tripple]

      std::size_t size() const { return map_.size(); }

      LocalKey& localKey ( std::size_t i ) { return map_[ i ]; }
      const LocalKey& localKey ( std::size_t i ) const { return map_[ i ]; }

    private:
      std::vector< LocalKey > map_;
    };

    struct BubbleElementDofMapperCodeFactory
    {
      // return the shape functions for a given reference element. If this
      // is not possible an empty DofMapperCode is returned.
      template< class Field, int dim >
      DofMapperCode operator() ( const ReferenceElement< Field, dim > &refElement ) const
      {
        if( refElement.type().isSimplex() )
          return compile( refElement, BubbleElementLocalKeyMap< dim >() );
        else
          return DofMapperCode();
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef SPACE_LOCALKEYMAP_HH
