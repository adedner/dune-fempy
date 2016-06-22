#include <config.h>

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

//#include "fullscheme.hh"
#include "uzawascheme.hh"
#include "prpscheme.hh"
#include "navierstokes.hh"
#include "testingquantities.hh"

// assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class HGridType>
double algorithm ( HGridType &grid, int step )
{
  // we want to solve the problem on the leaf elements of the grid
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > GridPartType;
  GridPartType gridPart(grid);

  //Function space functions map
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld, HGridType::dimensionworld+1 > FunctionSpaceType;// from R^2 -> R^3
  //typedef Dune::Fem::FunctionSpace< double, double,HGridType::dimensionworld, HGridType::dimensionworld > VelocityFunctionSpaceType;  //i.e velocity space functions mapping R^2 -> R^2
  //typedef Dune::Fem::FunctionSpace< double, double,HGridType::dimensionworld, 1 > PressureFunctionSpaceType;  //pressure space functions mapping R^2 -> R

  // create time provider and problem type
  Dune::Fem::GridTimeProvider< HGridType > timeProvider( grid );
  // typedef TemporalProblemInterface< FunctionSpaceType,VelocityFunctionSpaceType > NavierStokesProblemType;
  typedef NavierStokesProblemInterface< FunctionSpaceType> NavierStokesProblemType;
  NavierStokesProblemType* problemPtr = 0 ;

  const double stokesTimeStepTheta= Dune::Fem::Parameter::getValue< double >("navierstokes.stokestimestepfactor",0.3);//"theta" from the notes
  const double stokesTimeStepFactor=1/stokesTimeStepTheta; //1/theta
  const double burgersTimeStepTheta = 1-2*stokesTimeStepTheta;//"theta_b" from the notes
  const double burgersTimeStepFactor=1/burgersTimeStepTheta; //1/theta_b
  const double implicitFactor =  Dune::Fem::Parameter::getValue< double >("navierstokes.implicitfactor",0.5);// =1 fully implicit diffusion, =0 fully explicit diffusion, "alpha" in notes
  const double viscosity =  Dune::Fem::Parameter::getValue< double >("navierstokes.viscosity",0.3);// "eta" in the notes
  const double viscosityStokes = implicitFactor*viscosity;
  const double viscosityBurgers = implicitFactor*viscosity;

  const std::string NavierStokesProblemNames[] = { "channelflow","vorticityflow","movingplate","couetteflow", "karmanstreet" };
  enum { ChannelFlowProblem,VorticityFlowProblem,MovingPlateProblem,CouetteFlowProblem,KarmanStreetProblem};
  const int NavierStokesProblemNumber = Dune::Fem::Parameter::getEnum("navierstokes.problem",NavierStokesProblemNames,ChannelFlowProblem);
  switch (NavierStokesProblemNumber)
    {
    case ChannelFlowProblem: problemPtr = new ChannelFlow<FunctionSpaceType>(timeProvider,viscosity,stokesTimeStepFactor,burgersTimeStepFactor); break;
    case VorticityFlowProblem: problemPtr = new VorticityFlow<FunctionSpaceType>(timeProvider,viscosity,stokesTimeStepFactor,burgersTimeStepFactor); break;
    case MovingPlateProblem: problemPtr = new MovingPlate<FunctionSpaceType>(timeProvider,viscosity,stokesTimeStepFactor,burgersTimeStepFactor); break;
    case CouetteFlowProblem: problemPtr = new CouetteFlow<FunctionSpaceType>(timeProvider,viscosity,stokesTimeStepFactor,burgersTimeStepFactor); break;
    case KarmanStreetProblem: problemPtr = new KarmanVortexStreet<FunctionSpaceType>(timeProvider,viscosity,stokesTimeStepFactor,burgersTimeStepFactor); break;
    default: problemPtr = new ChannelFlow<FunctionSpaceType>(timeProvider,viscosity,stokesTimeStepFactor,burgersTimeStepFactor); break;
    }

  NavierStokesProblemType& problem = *problemPtr ;

  typedef Dune::Fem::GridFunctionAdapter< NavierStokesProblemType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problem, gridPart, 5 );


  // ------------------------------Stokes scheme initialisation
  typedef UzawaScheme< GridPartType > StokesSchemeType;
  StokesSchemeType StokesScheme( gridPart, problem, viscosityStokes,stokesTimeStepFactor );
  //--------------------------------Stokes finish initialisation

  //--------------------------------Burgers scheme initialisation
  typedef PRPScheme<GridPartType > BurgersSchemeType;
  BurgersSchemeType BurgersScheme( gridPart,problem,burgersTimeStepFactor,viscosityBurgers);//careful parameter order is different to in stokesscheme
  //--------------------------------Burgers finish initialistaion

  //--------------------------------Data output
  //typedef Dune::tuple< const typename StokesSchemeType::DiscreteFunctionType *,const typename BurgersSchemeType::DiscreteFunctionType *, const GridExactSolutionType * > IOTupleType;
  typedef Dune::tuple< const typename StokesSchemeType::DiscreteFunctionType *, const GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  //IOTupleType ioTuple( &(StokesScheme.solution()),&(BurgersScheme.solution()), &gridExactSolution); // tuple with pointers
  IOTupleType ioTuple( &(StokesScheme.solution()),&gridExactSolution); // tuple with pointers
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

 //--------------------------------time loop parameters
 const double endTime  = Dune::Fem::Parameter::getValue< double >( "navierstokes.endtime", 1.0 );
 const double dtreducefactor = Dune::Fem::Parameter::getValue< double >("navierstokes.reducetimestepfactor", 1 );
 double timeStep = Dune::Fem::Parameter::getValue< double >( "navierstokes.timestep", 0.125 );

 //typedef Dune::Fem::L2Norm< GridPartType > LTwoNormType;
 //LTwoNormType lTwoNorm( gridPart );
 double* Test;
 //----Test filewriting
  std::ofstream outfile;
  std::ostringstream s;
   s << "../navier_stokes_output/Test-" << step << ".txt";
  const std::string& tmp = s.str();
  outfile.open(tmp.c_str());

  outfile << "Time" << " "
      << "L^2 Error - Velocity"  << " "
      << "H^1 Error - Velocity" << " "
      << std::endl;

 timeStep *= pow(dtreducefactor,step);
 timeProvider.init( timeStep ) ;

 StokesScheme.initialize();
 Test = testintegral(gridExactSolution,StokesScheme.solution(),problem);
 dataOutput.write( timeProvider );

 //lTwoError =  lTwoNorm.distance( gridExactSolution, StokesScheme.solution() );
 std::cout << "L^2 Error to exact solution: " << *Test <<std::endl;
 // hOneError =  hOneNorm.distance( gridExactSolution, StokesScheme.solution() );
 //std::cout << "H^1 Error to exact solution: " << hOneError <<std::endl;
 outfile << timeProvider.time()<<"      "
     << *Test  <<"     " <<std::endl;
   //<< hOneError  <<std::endl;

  timeProvider.next( timeStep );
 //-------------------------------Begin solver
 for( ; timeProvider.time() < endTime; timeProvider.next( timeStep ) )
   {
     std::cout << "Prepare step 1 - Stokes" << std::endl;
     //Solve first step (prev step is stokes so no need to update)
     StokesScheme.preparestep1();
     std::cout << "Solve step 1 - Stokes" << std::endl;
     StokesScheme.solve( (timeProvider.time()<=timeStep) );

     //Solve second step (requires vel&pres update from stokes)
     std::cout << "Update Pressure" << std::endl;
     BurgersScheme.updatepressure(StokesScheme.pressure());
     std::cout << "Update Velocity" << std::endl;
     BurgersScheme.updatevelocity(StokesScheme.velocity());
     std::cout << "Prepare step 2 - Burgers" << std::endl;
     BurgersScheme.prepare();
     std::cout << "Solve step 2 - Burgers" << std::endl;
     BurgersScheme.solve((timeProvider.time()<=timeStep));

     //solve third step (requires vel update from burgers)
     std::cout << "Update Velocity" << std::endl;
     StokesScheme.updatevelocity(BurgersScheme.solution());
     std::cout << "Prepare step 3 - Stokes" << std::endl;
     StokesScheme.preparestep3();
     std::cout << "Solve step 3 - Stokes" << std::endl;
     StokesScheme.solve(false );
     dataOutput.write(timeProvider);

     //error calculation
     Test = testintegral(gridExactSolution,StokesScheme.solution(),problem);
     //lTwoError =  lTwoNorm.distance( gridExactSolution, StokesScheme.solution() );
     std::cout << "L^2 Error to exact solution: " << *Test <<std::endl;
     //  hOneError =  hOneNorm.distance( gridExactSolution, StokesScheme.solution() );
     //std::cout << "H^1 Error to exact solution: " << hOneError <<std::endl;
     outfile << timeProvider.time()<<"      "
         << *Test <<"     " <<std::endl;
       // << hOneError  <<std::endl;
   }

 outfile.close();
 return *Test ;
}

// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append overloaded parameters from the command line
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );

  // append default parameter file
  Dune::Fem::Parameter::append( "../navier_stokes_data/parameter" );

  // type of hierarchical grid
  typedef Dune::GridSelector::GridType  HGridType ;

  // create grid from DGF file
  const std::string gridkey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string gridfile = Dune::Fem::Parameter::getValue< std::string >( gridkey );

  // the method rank and size from MPIManager are static
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;

  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  // do initial load balance
  grid.loadBalance();

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "navierstokes.initialRefinements" );

  // number of global refinements to bisect grid width
  const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

  // refine grid
  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "navierstokes.repeats", 0 );

  // calculate first step
  double oldError = algorithm( grid, (repeats > 0) ? 0 : -1 );

  for( int step = 1; step <= repeats; ++step )
  {
    // refine globally such that grid with is bisected
    // and all memory is adjusted correctly
    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

    const double newError = algorithm( grid, step );
    const double eoc = log( oldError / newError ) / M_LN2;
    if( Dune::Fem::MPIManager::rank() == 0 )
    {
      std::cout << "Size: " << grid.size(0) << std::endl;
      std::cout << "Error: " << newError << std::endl;
      std::cout << "EOC( " << step << " ) = " << eoc << std::endl;
    }
    oldError = newError;
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
