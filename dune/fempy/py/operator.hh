#ifndef DUNE_FEMPY_PY_OPERATOR_HH
#define DUNE_FEMPY_PY_OPERATOR_HH

#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/function/discrete.hh>
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // registerOperatorJacobian
      // ------------------------

      template< class Operator, class... options >
      inline static auto registerOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
        -> void_t< typename Operator::JacobianRangeType >
      {
        typedef typename Operator::DomainFunctionType DomainFunction;

        using pybind11::operator""_a;

        cls.def( "jacobian", [] ( Operator &self, pybind11::object u, typename Operator::JacobianRangeType &jOp ) {
            asGridFunction< DomainFunction >( jOp.domainSpace().gridPart(), u, [ &self, &jOp ] ( const auto &u ) -> void_t< decltype( self.jacobian( u, jOp ) ) > { self.jacobian( u, jOp ); } );
          }, "u"_a, "jOp"_a );
      }

      template< class Operator, class... options >
      inline static void registerOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}

      template< class Operator, class... options >
      inline static void registerOperatorJacobian ( pybind11::class_< Operator, options... > cls )
      {
        registerOperatorJacobian ( cls, PriorityTag< 42 >() );
      }



      // registerOperator
      // ----------------

      template< class Operator, class... options >
      inline static void registerOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls )
      {
        typedef typename Operator::DomainFunctionType DomainFunction;

        using pybind11::operator""_a;

        cls.def( "__call__", [] ( Operator &self, pybind11::object u, RangeFunction &w ) {
            asGridFunction< DomainFunction >( w.gridPart(), u, [ &self, &w ] ( const auto &u ) -> void_t< decltype( self( u, w ) ) > { self( u, w ); } );
          }, "u"_a, "w"_a );

        registerOperatorJacobian( cls );
      }

    } // namespace detail



    // registerOperator
    // ----------------

    template< class Operator, class... options >
    inline static void registerOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls )
    {
      detail::registerOperator< Operator >( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_OPERATOR_HH
