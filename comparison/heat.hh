// Model
// -----

template< class GridPart, class Coefficient0 >
struct Model
{
  typedef GridPart GridPartType;
  typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPart::IntersectionType IntersectionType;
  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, GridPartType::dimensionworld, 1 > FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;
  static const int dimDomain = GridPartType::dimensionworld;
  static const int dimRange = 1;
  static const int dimLocal = GridPartType::dimension;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, std::tuple< std::shared_ptr< double >, std::shared_ptr< double >, std::shared_ptr< double > > >::element_type;
  static const std::size_t numConstants = 3;
  typedef std::tuple< Dune::Fem::FunctionSpace< DomainFieldType, double, dimDomain, 1 > > CoefficientFunctionSpaceTupleType;
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, std::tuple< Coefficient0 > >;

  Model ( const Coefficient0 &coefficient0 )
    : coefficients_( coefficient0 )
  {
    std::get< 0 >( constants_ ) = std::make_shared< double >();
    std::get< 1 >( constants_ ) = std::make_shared< double >();
    std::get< 2 >( constants_ ) = std::make_shared< double >();
  }

  bool init ( const EntityType &entity ) const
  {
    {
      entity_ = &entity;
      std::get< 0 >( coefficients_ ).init( entity );
    }
    return true;
  }

  const EntityType &entity () const
  {
    return *entity_;
  }

  std::string name () const
  {
    return "Model";
  }
  typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
  static const bool symmetric = false;

  template< class Point >
  void source ( const Point &x, const RangeType &u, const JacobianRangeType &du, RangeType &result ) const
  {
    using std::pow;
    using std::sin;
    const auto tmp0 = [ this ] ( auto x ) {
        typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::RangeType result;
        std::get< 0 >( coefficients_ ).evaluate( x, result );
        return result;
      }( x );
    const auto tmp1 = -1 * tmp0[ 0 ];
    const auto tmp2 = u[ 0 ] + tmp1;
    double tmp3 = constant< 2 >();
    double tmp4 = constant< 0 >();
    const auto tmp5 = -1 * tmp4;
    const auto tmp6 = 1.0 + tmp5;
    auto tmp7 = entity().geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp8 = tmp7[ 0 ] * tmp7[ 0 ];
    const auto tmp9 = tmp7[ 1 ] * tmp7[ 1 ];
    const auto tmp10 = tmp8 + tmp9;
    const auto tmp11 = 6.28318530718 * tmp10;
    double tmp12 = constant< 1 >();
    const auto tmp13 = std::pow( tmp12, 2 );
    const auto tmp14 = tmp11 * tmp13;
    const auto tmp15 = std::sin( tmp14 );
    const auto tmp16 = -1 * tmp15;
    const auto tmp17 = tmp6 * tmp16;
    const auto tmp18 = tmp3 + tmp12;
    const auto tmp19 = std::pow( tmp18, 2 );
    const auto tmp20 = tmp11 * tmp19;
    const auto tmp21 = std::sin( tmp20 );
    const auto tmp22 = -1 * tmp21;
    const auto tmp23 = tmp4 * tmp22;
    const auto tmp24 = tmp17 + tmp23;
    const auto tmp25 = tmp3 * tmp24;
    const auto tmp26 = tmp2 + tmp25;
    result[ 0 ] = tmp26;
  }

  template< class Point >
  void linSource ( const RangeType &ubar, const JacobianRangeType &dubar, const Point &x, const RangeType &u, const JacobianRangeType &du, RangeType &result ) const
  {
    result[ 0 ] = u[ 0 ];
  }

  template< class Point >
  void diffusiveFlux ( const Point &x, const RangeType &u, const JacobianRangeType &du, JacobianRangeType &result ) const
  {
    double tmp0 = constant< 2 >();
    double tmp1 = constant< 0 >();
    const auto tmp2 = (du[ 0 ])[ 0 ] * tmp1;
    const auto tmp3 = [ this ] ( auto x ) {
        typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::JacobianRangeType result;
        std::get< 0 >( coefficients_ ).jacobian( x, result );
        return result;
      }( x );
    const auto tmp4 = -1 * tmp1;
    const auto tmp5 = 1.0 + tmp4;
    const auto tmp6 = (tmp3[ 0 ])[ 0 ] * tmp5;
    const auto tmp7 = tmp2 + tmp6;
    const auto tmp8 = tmp0 * tmp7;
    const auto tmp9 = (du[ 0 ])[ 1 ] * tmp1;
    const auto tmp10 = (tmp3[ 0 ])[ 1 ] * tmp5;
    const auto tmp11 = tmp9 + tmp10;
    const auto tmp12 = tmp0 * tmp11;
    (result[ 0 ])[ 0 ] = tmp8;
    (result[ 0 ])[ 1 ] = tmp12;
  }

  template< class Point >
  void linDiffusiveFlux ( const RangeType &ubar, const JacobianRangeType &dubar, const Point &x, const RangeType &u, const JacobianRangeType &du, JacobianRangeType &result ) const
  {
    double tmp0 = constant< 2 >();
    double tmp1 = constant< 0 >();
    const auto tmp2 = (du[ 0 ])[ 0 ] * tmp1;
    const auto tmp3 = tmp0 * tmp2;
    const auto tmp4 = (du[ 0 ])[ 1 ] * tmp1;
    const auto tmp5 = tmp0 * tmp4;
    (result[ 0 ])[ 0 ] = tmp3;
    (result[ 0 ])[ 1 ] = tmp5;
  }

  template< class Point >
  void fluxDivergence ( const Point &x, const RangeType &u, const JacobianRangeType &du, const HessianRangeType &d2u, RangeType &result ) const
  {
    using std::pow;
    using std::sin;
    const auto tmp0 = [ this ] ( auto x ) {
        typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::RangeType result;
        std::get< 0 >( coefficients_ ).evaluate( x, result );
        return result;
      }( x );
    const auto tmp1 = -1 * tmp0[ 0 ];
    const auto tmp2 = u[ 0 ] + tmp1;
    double tmp3 = constant< 2 >();
    double tmp4 = constant< 0 >();
    const auto tmp5 = -1 * tmp4;
    const auto tmp6 = 1.0 + tmp5;
    auto tmp7 = entity().geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp8 = tmp7[ 0 ] * tmp7[ 0 ];
    const auto tmp9 = tmp7[ 1 ] * tmp7[ 1 ];
    const auto tmp10 = tmp8 + tmp9;
    const auto tmp11 = 6.28318530718 * tmp10;
    double tmp12 = constant< 1 >();
    const auto tmp13 = std::pow( tmp12, 2 );
    const auto tmp14 = tmp11 * tmp13;
    const auto tmp15 = std::sin( tmp14 );
    const auto tmp16 = -1 * tmp15;
    const auto tmp17 = tmp6 * tmp16;
    const auto tmp18 = tmp3 + tmp12;
    const auto tmp19 = std::pow( tmp18, 2 );
    const auto tmp20 = tmp11 * tmp19;
    const auto tmp21 = std::sin( tmp20 );
    const auto tmp22 = -1 * tmp21;
    const auto tmp23 = tmp4 * tmp22;
    const auto tmp24 = tmp17 + tmp23;
    const auto tmp25 = tmp3 * tmp24;
    const auto tmp26 = tmp2 + tmp25;
    const auto tmp27 = ((d2u[ 0 ])[ 0 ])[ 0 ] * tmp4;
    const auto tmp28 = [ this ] ( auto x ) {
        typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::HessianRangeType result;
        std::get< 0 >( coefficients_ ).hessian( x, result );
        return result;
      }( x );
    const auto tmp29 = ((tmp28[ 0 ])[ 0 ])[ 0 ] * tmp6;
    const auto tmp30 = tmp27 + tmp29;
    const auto tmp31 = tmp3 * tmp30;
    const auto tmp32 = ((d2u[ 0 ])[ 1 ])[ 1 ] * tmp4;
    const auto tmp33 = ((tmp28[ 0 ])[ 1 ])[ 1 ] * tmp6;
    const auto tmp34 = tmp32 + tmp33;
    const auto tmp35 = tmp3 * tmp34;
    const auto tmp36 = tmp31 + tmp35;
    const auto tmp37 = -1 * tmp36;
    const auto tmp38 = tmp26 + tmp37;
    result[ 0 ] = tmp38;
  }

  template< class Point >
  void alpha ( const Point &x, const RangeType &u, RangeType &result ) const
  {
    result[ 0 ] = 0;
  }

  template< class Point >
  void linAlpha ( const RangeType &ubar, const Point &x, const RangeType &u, RangeType &result ) const
  {
    result[ 0 ] = 0;
  }

  bool hasDirichletBoundary () const
  {
    return false;
  }

  bool hasNeumanBoundary () const
  {
    return false;
  }

  bool isDirichletIntersection ( const IntersectionType &intersection, Dune::FieldVector< int, dimRange > &dirichletComponent ) const
  {
    return false;
  }

  template< class Point >
  void dirichlet ( int bndId, const Point &x, RangeType &result ) const
  {
    result = RangeType( 0 );
  }

  template< std::size_t i >
  const ConstantType< i > &constant () const
  {
    return *std::get< i >( constants_ );
  }

  template< std::size_t i >
  ConstantType< i > &constant ()
  {
    return *std::get< i >( constants_ );
  }

  template< std::size_t i >
  const CoefficientType< i > &coefficient () const
  {
    return std::get< i >( coefficients_ );
  }

  template< std::size_t i >
  CoefficientType< i > &coefficient ()
  {
    return std::get< i >( coefficients_ );
  }
  class BoundaryWrapper
  {
    const Model& impl_;
    int bndId_;
    public:
    BoundaryWrapper( const Model& impl, int bndId )
    : impl_( impl ), bndId_(bndId) {}

    //! evaluate function
    template <class Point>
    void evaluate( const Point& x, RangeType& ret ) const
    {
      impl_.dirichlet(bndId_,Dune::Fem::coordinate(x),ret);
    }
    //! jacobian function (only for exact)
    void jacobian( const DomainType& x, JacobianRangeType& ret ) const
    {
      DUNE_THROW(Dune::NotImplemented,"rhs jacobian not implemented");
    }
  };

private:
  mutable const EntityType *entity_ = nullptr;
  mutable std::tuple< std::shared_ptr< double >, std::shared_ptr< double >, std::shared_ptr< double > > constants_;
  mutable std::tuple< Coefficient0 > coefficients_;
};
