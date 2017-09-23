#include <config.h>

// iostream includes
#include <iostream>

#include <dune/fem/schemes/elliptic.hh>
#include <dune/fem/schemes/femscheme.hh>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/io/file/dataoutput.hh>

// move to fem?
#include <dune/fempy/grid/gridpartadapter.hh>

// include header for heat model (generated from ufl)
#include "heat.hh"

const double dt = 0.02;

template <class GridPart>
struct Initial : public Dune::Fem::BindableGridFunction< GridPart, Dune::Dim<1> >
{
  typedef Dune::Fem::BindableGridFunction<GridPart, Dune::Dim<1> > Base;
  using Base::Base;

  template <class Point>
  void evaluate(const Point &xhat, typename Base::RangeType &ret) const
  {
    auto x = Base::global(xhat);
    ret[0] = std::atan( pow(10.0 * x[0] * (1-x[0]) * x[1] * (1-x[1]),2) );
  }
  unsigned int order() const { return 5; }
  std::string name() const { return "Initial"; }
};


template <class HGridType>
void algorithm ( HGridType &grid, int step )
{
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld, 1 > FunctionSpaceType;
  typedef Dune::FemPy::GridPartAdapter<typename HGridType::LeafGridView> GridPartType;
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 2> SpaceType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< SpaceType > DiscreteFunctionType;
  typedef Dune::Fem::ConstLocalDiscreteFunction<DiscreteFunctionType> LocalFunctionType;
  typedef Model< GridPartType, LocalFunctionType > ModelType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef DifferentiableEllipticOperator< LinearOperatorType, ModelType > OperatorType;
  typedef Dune::Fem::KrylovInverseOperator< DiscreteFunctionType > InverseLinearOperatorType;
  typedef FemScheme< OperatorType, InverseLinearOperatorType> SchemeType;

  auto gridView = grid.leafGridView();
  GridPartType gridPart(gridView);
  SpaceType space( gridPart );
  DiscreteFunctionType solution( "solution", space );
  Dune::Fem::interpolate(Initial<GridPartType>(gridPart),solution);
  DiscreteFunctionType u_n(  "previous", space );
  LocalFunctionType lu_n(u_n);
  ModelType model( lu_n );
  // t=1, theta=0, dt=2
  model.template constant<2>() = dt;
  model.template constant<0>() = 0.5;
  SchemeType scheme( space, model );

  typedef std::tuple< DiscreteFunctionType * > IOTupleType;
  IOTupleType ioTuple( &solution );
  Dune::Fem::DataOutput<HGridType,IOTupleType> dataOutput( grid, ioTuple );

  dataOutput.writeData( 0 );
  const int steps = int(1. / dt);
  for (int n=0;n<steps;++n)
  {
    model.template constant<1>() = n*dt;
    u_n.assign(solution);
    scheme.solve(solution);
  }
  dataOutput.writeData( steps );
}

// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );
  Dune::Fem::Parameter::append( "parameter" );

  // type of hierarchical grid
  typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > HGridType;

  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( "unitcube.dgf" );
  HGridType& grid = *gridPtr ;

  // do initial load balance
  grid.loadBalance();

  algorithm( grid, 0 );

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
