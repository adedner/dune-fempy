#ifndef DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
#define DUNE_FEMPY_PY_DISCRETEFUNCTION_HH

#include <dune/fem/space/common/interpolate.hh>

#include <dune/corepy/pybind11/pybind11.h>
#include <dune/corepy/pybind11/extensions.h>

#include <dune/fempy/py/function/grid.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/py/grid/restrictprolong.hh>

#include <dune/fem/function/vectorfunction/vectorfunction.hh>
#include <dune/fempy/py/common/numpyvector.hh>

#include <dune/corepy/pybind11/eigen.h>

namespace Dune
{

  namespace FemPy
  {

    // interpolant
    // -----------
    template< class DF, class GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value, int > = 0 >
    void interpolant ( DF &df, std::string name, const GF &gf )
    {
      Fem::interpolate( gf, df );
    }

    template< class DF >
    void interpolant ( DF &df, typename DF::RangeType value )
    {
      return interpolant( df, simpleGridFunction( df.space().gridPart(), [ value ] ( typename DF::DomainType ) { return value; }, 0 ) );
    }


    // registerDiscreteFunction
    // -------------

    namespace detail
    {
      template< class DF, class Cls >
      void registerDFConstructor( Cls &cls, std::false_type )
      {}
      template<class DF, class Cls >
      void registerDFConstructor( Cls &cls, std::true_type )
      {
        typedef typename DF::DiscreteFunctionSpaceType Space;
        cls.def( "__init__", [] ( DF &instance, Space &space, const std::string &name ) {
            new( &instance ) DF( std::move(name), space );
          }, pybind11::keep_alive< 1, 2 >() );
      }

#if 0
      template <class DF>
      auto registerDFBuffer(DF &instance,int)
      -> decltype(instance.dofVector().array().data(),pybind11::buffer_info())
      {
        typedef typename DF::RangeType Value;
        typedef typename Value::field_type FieldType;
        return pybind11::buffer_info(
            instance.dofVector().array().data(),    /* Pointer to buffer */
            sizeof(FieldType),                      /* Size of one scalar */
            pybind11::format_descriptor<FieldType>::format(), /* Python struct-style format descriptor */
            1,                                      /* Number of dimensions */
            { instance.size() },                    /* Buffer dimensions */
            { sizeof(FieldType) }                   /* Strides (in bytes) for each index */
        );
      }
      template <class DF>
      pybind11::buffer_info registerDFBuffer(DF &instance,long)
      {
        // THROW SOME EXCEPTION
        typedef typename DF::RangeType Value;
        typedef typename Value::field_type FieldType;
        return pybind11::buffer_info(
            nullptr,                                /* Pointer to buffer */
            sizeof(FieldType),                      /* Size of one scalar */
            pybind11::format_descriptor<FieldType>::format(), /* Python struct-style format descriptor */
            1,                                      /* Number of dimensions */
            { instance.size() },                    /* Buffer dimensions */
            { sizeof(FieldType) }                   /* Strides (in bytes) for each index */
        );
      }
#endif
      template <class DF, class Cls>
      auto registerDFBuffer(Cls &cls,int)
      -> decltype(std::declval<DF>().dofVector().array().data(),void())
      {
        typedef typename DF::RangeType Value;
        typedef typename Value::field_type FieldType;
        cls.def_buffer( [](DF &instance) -> pybind11::buffer_info {
            return pybind11::buffer_info(
                instance.dofVector().array().data(),    /* Pointer to buffer */
                sizeof(FieldType),                      /* Size of one scalar */
                pybind11::format_descriptor<FieldType>::format(), /* Python struct-style format descriptor */
                1,                                      /* Number of dimensions */
                { instance.size() },                    /* Buffer dimensions */
                { sizeof(FieldType) }                   /* Strides (in bytes) for each index */
            );
          }); // ????  pybind11::keep_alive<0,1>() );
      }
      template <class DF,class Cls>
      void registerDFBuffer(Cls &cls,long)
      {
      }

      template< class DF, class Cls >
      void registerDiscreteFunction ( pybind11::module module, Cls &cls )
      {
        typedef typename DF::DiscreteFunctionSpaceType Space;
        typedef typename DF::GridPartType GridPart;
        typedef typename DF::RangeType Value;
        typedef typename GridPart::GridType Grid;
        typedef typename Value::field_type FieldType;

        detail::registerGridFunction< DF >( module, cls);

        detail::clsVirtualizedGridFunction< GridPart, Value >( module ).def( "__init__", [] ( VirtualizedGridFunction< GridPart, Value > &instance, DF &df ) {
            new (&instance) VirtualizedGridFunction< GridPart, Value >( pyGridFunction( df ) );
          } );
        pybind11::implicitly_convertible< DF, VirtualizedGridFunction< GridPart, Value > >();

        detail::clsVirtualizedRestrictProlong< Grid >( module ).def( "__init__", [] ( VirtualizedRestrictProlong< Grid > &instance, DF &df ) {
            new (&instance) VirtualizedRestrictProlong< Grid >( df );
          }, pybind11::keep_alive< 1, 2 >() );
        pybind11::implicitly_convertible< DF, VirtualizedRestrictProlong< Grid > >();

        cls.def_property_readonly( "space", [](DF &df) -> const typename DF::DiscreteFunctionSpaceType& {return df.space();} );
        cls.def("clear", [] (DF &instance) { instance.clear(); } );

        registerDFConstructor< DF, Cls >( cls, std::is_constructible< DF, const std::string&, const Space& >() );

        cls.def( "assign", [] ( DF &instance, const DF &other ) { instance.assign(other); } );

        typedef VirtualizedGridFunction< GridPart, typename Space::RangeType > GridFunction;
        cls.def( "interpolate", [] ( DF &df, const GridFunction &gf ) {
            Fem::interpolate( gf, df );
          } );
        cls.def( "interpolate", [] ( DF &df, typename Space::RangeType value ) {
            Fem::interpolate(
                  simpleGridFunction( df.space().gridPart(), [ value ] ( typename DF::DomainType ) { return value; },0
            ), df );
          } );

        typedef typename DF::DofVectorType DofVector;
        if (!pybind11::already_registered<DofVector>())
        {
          auto clsDof = pybind11::class_<DofVector>( module, "DofVector");
        }
        cls.def( "dofVector", [] ( DF &instance ) -> DofVector&{ return instance.dofVector(); },
                 pybind11::return_value_policy::reference_internal );
        cls.def( "assign", [] ( DF &instance, const DofVector &other ) { instance.dofVector() = other; } );

        registerDFBuffer<DF>(cls,0);
      }
    }
    template< class DF, class Holder, class AliasType >
    void registerDiscreteFunction ( pybind11::module module, pybind11::class_<DF,Holder,AliasType> &cls )
    {
      detail::registerDiscreteFunction<DF>(module,cls);
    }
    template< class Space, class Field, class Holder, class AliasType >
    void registerDiscreteFunction ( pybind11::module module,
          pybind11::class_<Dune::Fem::VectorDiscreteFunction<Space,Dune::FemPy::NumPyVector<Field>>,
          Holder,AliasType> &cls )
    {
      typedef Dune::Fem::VectorDiscreteFunction<Space,Dune::FemPy::NumPyVector<Field>> DF;
      detail::registerDiscreteFunction<DF>(module,cls);
      typedef typename DF::VectorType VectorType;
      cls.def( "__init__", [] ( DF &instance, Space &space, const std::string &name, pybind11::buffer &dof ) {
          VectorType *vec = new VectorType(dof);
          new( &instance ) DF( std::move(name), space, *vec );
        }, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 4 >() );
    }
  } // namespace FemPy

} // namespace Dune


#endif // DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
