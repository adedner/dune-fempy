#ifndef DUNE_FEMPY_COMMON_UTILITY_HH
#define DUNE_FEMPY_COMMON_UTILITY_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/std/type_traits.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typeutilities.hh>

namespace Dune
{

  namespace FemPy
  {

    // same_type
    // ---------

    template< class T, class... U >
    using same_type = std::enable_if_t< Std::conjunction< std::is_same< T, U >... >::value, T >;



    // FirstMatch
    // ----------

    template< class... Fs >
    class FirstMatch
    {
      template< class... Args >
      struct CallTraits
      {
        template< class F >
        static decltype( std::declval< F >()( std::declval< Args >()... ), std::true_type() ) __Matches ( PriorityTag< 1 > );

        template< class F >
        static std::false_type __Matches ( PriorityTag< 0 > );

      public:
        template< class F >
        using Matches = decltype( __Matches< F >( PriorityTag< 42 >() ) );
      };

      template< class... Args >
      using ConstFirstMatchIndex = FirstPredicateIndex< std::tuple< const Fs &... >, CallTraits< Args... >::template Matches >;

      template< class... Args >
      using FirstMatchIndex = FirstPredicateIndex< std::tuple< Fs &... >, CallTraits< Args... >::template Matches >;

    public:
      explicit FirstMatch ( Fs... fs ) : fs_( std::move( fs )... ) {}

      template< class... Args >
      decltype( auto ) operator() ( Args &&... args ) const
      {
        std::get< ConstFirstMatchIndex< Args... >::value >( fs_ )( std::forward< Args >( args )... );
      }

      template< class... Args >
      decltype( auto ) operator() ( Args &&... args )
      {
        std::get< FirstMatchIndex< Args... >::value >( fs_ )( std::forward< Args >( args )... );
      }

    private:
      std::tuple< Fs... > fs_;
    };



    // firstMatch
    // ----------

    template< class... Fs >
    inline static FirstMatch< std::decay_t< Fs >... > firstMatch ( Fs &&... fs )
    {
      return FirstMatch< std::decay_t< Fs >... >( std::forward< Fs >( fs )... );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_COMMON_UTILITY_HH
