#ifndef DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH
#define DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH

#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/operator/common/temporarylocalmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalMatrixEntry
    // ----------------

    template< class LocalMatrix >
    class LocalMatrixEntry
    {
      typedef LocalMatrixEntry< LocalMatrix > ThisType;

    public:
      typedef LocalMatrix LocalMatrixType;

      typedef typename LocalMatrixType::RangeFieldType RangeFieldType;

      LocalMatrixEntry ( LocalMatrixType &localMatrix, unsigned int row, unsigned int col )
        : localMatrix_( localMatrix ), row_( row ), col_( col )
      {}

      operator RangeFieldType () const { return localMatrix_.get( row_, col_ ); }

      ThisType &operator= ( const RangeFieldType &value ) { localMatrix_.set( row_, col_, value ); return *this; }

      ThisType &operator+= ( const RangeFieldType &value ) { localMatrix_.add( row_, col_, value ); return *this; }
      ThisType &operator-= ( const RangeFieldType &value ) { localMatrix_.add( row_, col_, -value ); return *this; }

      ThisType &operator*= ( const RangeFieldType &value ) { localMatrix_.set( row_, col_, localMatrix_.get( row_, col_ ) * value ); return *this; }
      ThisType &operator/= ( const RangeFieldType &value ) { localMatrix_.set( row_, col_, localMatrix_.get( row_, col_ ) / value ); return *this; }

    private:
      LocalMatrixType &localMatrix_;
      unsigned int row_, col_;
    };



    // LocalMatrixColumn
    // -----------------

    template< class LocalMatrix, class = void >
    struct LocalMatrixColumn
    {
      typedef LocalMatrix LocalMatrixType;

      typedef typename LocalMatrixType::RangeFieldType RangeFieldType;
      typedef typename LocalMatrixType::RangeBasisFunctionSetType BasisFunctionSetType;

      typedef typename BasisFunctionSetType::RangeType RangeType;
      typedef typename BasisFunctionSetType::JacobianRangeType JacobianRangeType;
      typedef typename BasisFunctionSetType::HessianRangeType HessianRangeType;

      typedef RangeFieldType value_type;
      typedef unsigned int size_type;

      LocalMatrixColumn ( LocalMatrixType &localMatrix, unsigned int col )
        : localMatrix_( localMatrix ), col_( col )
      {}

      RangeFieldType operator[] ( size_type row ) const { return localMatrix_.get( row, col_ ); }

      LocalMatrixEntry< LocalMatrixType > operator[] ( size_type row ) { return LocalMatrixEntry< LocalMatrixType >( localMatrix_, row, col_ ); }

      template< class Point, class... Factor >
      void axpy ( const Point &x, Factor &&... factor )
      {
        basisFunctionSet().axpy( x, std::forward< Factor >( factor )..., *this );
      }

      const BasisFunctionSetType &basisFunctionSet () const { return localMatrix_.rangeBasisFunctionSet(); }

    private:
      LocalMatrixType &localMatrix_;
      unsigned int col_;
    };



    // LocalMatrixColumn for Indexable LocalMatrix
    // -------------------------------------------

    template< class LocalMatrix >
    struct LocalMatrixColumn< LocalMatrix, std::enable_if_t< is_indexable< LocalMatrix, unsigned int >::value > >
    {
      typedef LocalMatrix LocalMatrixType;

      typedef typename LocalMatrixType::RangeBasisFunctionSetType BasisFunctionSetType;

      typedef typename BasisFunctionSetType::RangeType RangeType;
      typedef typename BasisFunctionSetType::JacobianRangeType JacobianRangeType;
      typedef typename BasisFunctionSetType::HessianRangeType HessianRangeType;

      typedef typename FieldTraits< RangeType >::field_type RangeFieldType;

      typedef std::decay_t< decltype( std::declval< const LocalMatrix & >()[ 0u ][ 0u ] ) > value_type;
      typedef unsigned int size_type;

      LocalMatrixColumn ( LocalMatrixType &localMatrix, unsigned int col )
        : localMatrix_( localMatrix ), col_( col )
      {}

      const value_type &operator[] ( size_type row ) const { return localMatrix_[ row ][ col_ ]; }
      value_type &operator[] ( size_type row ) { return localMatrix_[ row ][ col_ ]; }

      template< class Point, class... Factor >
      void axpy ( const Point &x, Factor &&... factor )
      {
        basisFunctionSet().axpy( x, std::forward< Factor >( factor )..., *this );
      }

      const BasisFunctionSetType &basisFunctionSet () const { return localMatrix_.rangeBasisFunctionSet(); }

    private:
      LocalMatrixType &localMatrix_;
      unsigned int col_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH
