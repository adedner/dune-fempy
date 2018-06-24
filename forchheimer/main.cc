#include <config.h>

// iostream includes
#include <iostream>
#include <complex>

#include <dune/grid/yaspgrid.hh>

#include <dune/fempy/grid/gridpartadapter.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/io/file/dataoutput.hh>

// include header of elliptic solver
#include <dune/fem/schemes/elliptic.hh>
#include <dune/fem/schemes/femscheme.hh>

// include generated model
#include <forchheimer/forchheimer.hh>

template <class GridPart>
struct Initial : public Dune::Fem::BindableGridFunction< GridPart, Dune::Dim<1> >
{
  typedef Dune::Fem::BindableGridFunction<GridPart, Dune::Dim<1> > Base;
  using Base::Base;
  template <class Point>
  void evaluate(const Point &xhat, typename Base::RangeType &ret) const
  {
    auto x = Base::global(xhat);
    ret[0] = 1./2.*x.two_norm2() - 1./3.*(pow(x[0],3) - pow(x[1],3)) + 1.;
  }
  unsigned int order() const { return 5; }
  std::string name() const { return "Initial"; }
};

int main ( int argc, char **argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );
  Dune::Fem::Parameter::append( argc, argv );
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );
  Dune::Fem::Parameter::append( "parameter" );
  typedef Dune::YaspGrid<2> HGridType ;

  const std::string gridkey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string gridfile = Dune::Fem::Parameter::getValue< std::string >( gridkey );
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  auto gridView = grid.leafGridView();
  Dune::FemPy::GridPartAdapter<decltype(gridView)> gridPart(gridView);
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld, 1 > FunctionSpaceType;
  Dune::Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType,decltype(gridPart),2> space(gridPart);
  Dune::Fem::AdaptiveDiscreteFunction<decltype(space)> solution("solution",space);
  decltype(solution) previous(solution);

  Dune::Fem::interpolate(Initial<decltype(gridPart)>(gridPart),solution);

  Model<decltype(gridPart),typename decltype(previous)::LocalFunctionType> model( previous.localFunction() );

  typedef FemScheme< DifferentiableEllipticOperator<
          Dune::Fem::SparseRowLinearOperator<decltype(solution),decltype(solution)>,decltype(model)>,
          Dune::Fem::KrylovInverseOperator<decltype(solution)> > SchemeType;
  SchemeType scheme( space, model );

  std::tuple< decltype(solution)* > ioTuple( &solution );
  Dune::Fem::DataOutput<HGridType,decltype(ioTuple)> dataOutput( grid, ioTuple );
  dataOutput.writeData( 0 );

  Dune::Fem::GridTimeProvider< HGridType > timeProvider( grid );
  double timeStep = 0.05;
  model.template constant<0>() = timeStep;    // model.dt
  for( timeProvider.init( timeStep ); timeProvider.time() < 1.0; timeProvider.next( timeStep ) )
  {
    previous.assign(solution);
    model.template constant<1>() = timeProvider.time();    // model.t
    scheme.solve( solution );
  }
  dataOutput.writeData( 1 );

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
