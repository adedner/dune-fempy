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

      // GeneralGridFunction
      // -------------------

      template< class DiscreteFunction >
      using GeneralGridFunction = VirtualizedGridFunction< typename DiscreteFunction::GridPartType, typename DiscreteFunction::RangeType >;



      // registerOperatorCall
      // --------------------

      template< class Operator, class... options >
      inline static auto registerOperatorCall ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< const Operator & >()( std::declval< const GeneralGridFunction< typename Operator::DomainFunctionType > & >(), std::declval< typename Operator::RangeFunctionType & >() ) >
      {
        typedef typename Operator::DomainFunctionType DomainFunction;

        using pybind11::operator""_a;

        cls.def( "__call__", [] ( Operator &self, pybind11::object u, RangeFunction &w ) {
            if( pybind11::isinstance< DomainFunction >( u ) )
              self( pybind11::cast< const DomainFunction & >( u ), w );
            else
              asGridFunction< typename DomainFunction::RangeType >( w.gridPart(), u, [ &self, &w ] ( const auto &u ) { self( u, w ); } );
          }, "u"_a, "w"_a );
      }

      template< class Operator, class... options >
      inline static void registerOperatorCall ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {
        typedef typename Operator::DomainFunctionType DomainFunction;

        using pybind11::operator""_a;

        cls.def( "__call__", [] ( Operator &self, const DomainFunction &u, RangeFunction &w ) { self( u, w ); }, "u"_a, "w"_a );
      }

      template< class Operator, class... options >
      inline static void registerOperatorCall ( pybind11::class_< Operator, options... > cls )
      {
        registerOperatorCall( cls, PriorityTag< 42 >() );
      }



      // registerOperatorJacobian
      // ------------------------

      template< class Operator, class... options >
      inline static void registerGeneralOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 2 > )
        -> void_t< decltype( std::declval< const Operator & >().jacobian( std::declval< const GeneralGridFunction< typename Operator::DomainFunctionType > & >(), std::declval< typename Operator::JacobianOperatorType >() ) >
      {
        typedef typename Operator::DomainFunctionType DomainFunction;

        using pybind11::operator""_a;

        cls.def( "jacobian", [] ( Operator &self, pybind11::object u, typename Operator::JacobianRangeType &jOp ) {
            if( pybind11::isinstance< DomainFunction >( u ) )
              self.jacobian( pybind11::cast< const DomainFunction & >( u ), jOp );
            else
              asGridFunction< typename DomainFunction::RangeType >( jOp.domainSpace().gridPart(), u, [ &self, &jOp ] ( const auto &u ) { self.jacobian( u, jOp ); } );
          }, "u"_a, "jOp"_a );
      }

      template< class Operator, class... options >
      inline static auto registerOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
        -> void_t< std::declval< const Operator & >().jacobian( std::declval< const typename Operator::DomainFunctionType & >(), std::declval< typename Operator::JacobianOperatorType >() ) >
      {
        typedef typename Operator::DomainFunctionType DomainFunction;

        using pybind11::operator""_a;

        cls.def( "jacobian", [] ( Operator &self, const DomainFunction &u, typename Operator::JacobianRangeType &jOp ) {
            self.jacobian( u, jOp );
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
        registerOperatorCall( cls );
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
