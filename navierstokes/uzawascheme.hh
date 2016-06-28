#ifndef STOKES_UZAWASCHEME_HH
#define STOKES_UZAWASCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/combinedspace.hh>

// adaptation ...
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/adaptmanager.hh>

// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>

// include linear operators
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>

#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/pardginverseoperators.hh>
#include <dune/fem/solver/oemsolver.hh>

// lagrange interpolation
#include <dune/fem/operator/lagrangeinterpolation.hh>

/*********************************************************/

// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>

// local includes
//#include "probleminterface.hh"
#include "navierstokes.hh"

#include "stokesmodel.hh"

#include "rhs.hh"
#include "elliptic.hh"
#include "noslipconstraints.hh"

#include "space.hh"

// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int step )
  : step_( step )
  {}

  DataOutputParameters ( const DataOutputParameters &other )
  : step_( other.step_ )
  {}

  std::string prefix () const
  {
    std::stringstream s;
    s << "navier_stokes-" << step_ << "-";
    return s.str();
  }
  std::string path () const
  {
    std::stringstream s;
    s << "../navier_stokes_output/";
    return s.str();
  }

private:
  int step_;
};

template < class VelocitySpace, class PressureSpace >
class UzawaScheme
{
public:
  typedef typename VelocitySpace::GridPartType GridPartType;

  typedef StokesMainModel<GridPartType>       MainModelType;
  typedef StokesGradModel<GridPartType>       GradModelType;
  typedef StokesDivergenceModel<GridPartType> DivergenceModelType;
  typedef StokesMassModel<GridPartType>       MassModelType;
  typedef StokesPrecondModel<GridPartType>    PrecondModelType;
  typedef StokesTransportModel<GridPartType>  TransportModelType;

  typedef typename GridPartType::GridType GridType;

  typedef typename GradModelType::VelocityFunctionSpaceType VelocityFunctionSpaceType;
  typedef typename GradModelType::PressureFunctionSpaceType PressureFunctionSpaceType;

#if MINIELEMENT
#error SHOULD NOT BE USING MINIELEMENT
 typedef Dune::Fem::BubbleElementSpace< VelocityFunctionSpaceType, GridPartType > VelocitySpaceType;
#else
 typedef VelocitySpace VelocitySpaceType;
#endif
  typedef PressureSpace PressureSpaceType;
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< VelocitySpaceType > VelocityDiscreteFunctionType;
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< PressureSpaceType > PressureDiscreteFunctionType;

  typedef Dune::Fem::ISTLLinearOperator< VelocityDiscreteFunctionType, VelocityDiscreteFunctionType > MainLinearOperatorType;
  typedef Dune::Fem::ISTLCGOp< VelocityDiscreteFunctionType, MainLinearOperatorType >              MainLinearInverseOperatorType;
  typedef Dune::Fem::ISTLLinearOperator< PressureDiscreteFunctionType, VelocityDiscreteFunctionType > GradLinearOperatorType;
  typedef Dune::Fem::ISTLLinearOperator< VelocityDiscreteFunctionType, PressureDiscreteFunctionType > DivergenceLinearOperatorType;
  typedef Dune::Fem::ISTLLinearOperator< PressureDiscreteFunctionType, PressureDiscreteFunctionType > MassLinearOperatorType;
  typedef Dune::Fem::ISTLCGOp< PressureDiscreteFunctionType, MassLinearOperatorType >                 MassLinearInverseOperatorType;
  typedef Dune::Fem::ISTLLinearOperator< PressureDiscreteFunctionType, PressureDiscreteFunctionType > PrecondLinearOperatorType;
  typedef Dune::Fem::ISTLCGOp< PressureDiscreteFunctionType, PrecondLinearOperatorType >           PrecondLinearInverseOperatorType;

  typedef EllipticOperator<  VelocityDiscreteFunctionType, VelocityDiscreteFunctionType,TransportModelType, NoConstraints > TransportOperatorType;

  typedef typename Dune::Fem::FunctionSpace<double,double,GridPartType::dimensionworld,GridPartType::dimensionworld+1> FullFunctionSpaceType;
  typedef NavierStokesProblemInterface< FullFunctionSpaceType> ProblemType ;

  typedef Dune::DirichletConstraints<MainModelType, VelocitySpaceType> MainConstraintsType;
  typedef DifferentiableEllipticOperator< MainLinearOperatorType, MainModelType, MainConstraintsType > MainOperatorType;
  typedef DifferentiableEllipticOperator< GradLinearOperatorType, GradModelType, NoConstraints > GradOperatorType;
  typedef DifferentiableEllipticOperator< DivergenceLinearOperatorType, DivergenceModelType, NoConstraints > DivergenceOperatorType;
  typedef DifferentiableEllipticOperator< MassLinearOperatorType, MassModelType, NoConstraints > MassOperatorType;
  typedef DifferentiableEllipticOperator< PrecondLinearOperatorType, PrecondModelType, NoConstraints > PrecondOperatorType;

  typedef typename  MainModelType::InitialVelocityFunctionType InitialVelocityFunctionType;
  typedef typename  MainModelType::InitialPressureFunctionType InitialPressureFunctionType;

  UzawaScheme( const VelocitySpace &velocitySpace, const PressureSpace &pressureSpace, const ProblemType &problem, double mu, double nu )
    : gridPart_( velocitySpace.gridPart() ),
      mu_(mu), nu_(nu),
      mainModel_( problem, gridPart_, mu, nu, true ),
      explicitMainModel_( problem, gridPart_, mu, nu, false ),
      gradModel_( gridPart_ ),
      divModel_( gridPart_ ),
      massModel_( gridPart_ ),
      precondModel_( problem, gridPart_, mu, nu ),
      transportModel_( problem, gridPart_ ),
      velocitySpace_( velocitySpace ),
      pressureSpace_( pressureSpace ),
      rhsU_( "rhsU", velocitySpace_ ),
      xi_( "xi", velocitySpace_ ),
      rhsP_( "rhsP", pressureSpace_ ),
      d_( "d", pressureSpace_ ),
      r_( "r", pressureSpace_ ),
      precond_( "precond", pressureSpace_ ),
      mainOperator_( mainModel_, velocitySpace_ ),
      explicitMainOperator_( explicitMainModel_, velocitySpace_ ),
      gradOperator_( gradModel_, velocitySpace_ ),
      divOperator_( divModel_, pressureSpace_ ),
      massOperator_( massModel_, pressureSpace_ ),
      precondOperator_( precondModel_, pressureSpace_ ),
      transportOperator_(transportModel_,velocitySpace_),
      // create linear operator (domainSpace,rangeSpace)
      mainLinearOperator_( "assembled main velocity operator", velocitySpace_, velocitySpace_ ),
      explicitMainLinearOperator_( "assembled main velocity operator", velocitySpace_, velocitySpace_ ),
      gradLinearOperator_( "assembled gradient operator", pressureSpace_, velocitySpace_ ),
      divLinearOperator_( "assembled divergence operator", velocitySpace_, pressureSpace_ ),
      massLinearOperator_( "assembled mass operator", pressureSpace_, pressureSpace_ ),
      precondLinearOperator_( "assembled precond operator", pressureSpace_, pressureSpace_ ),
      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "stokes.solvereps", 1e-8 ) ),
      usePrecond_( Dune::Fem::Parameter::getValue< bool >( "stokes.cg.usepreconditioner", true ) )
  {
    // set all DoF to zero
    //pressure_.clear();
    //velocity_.clear();
  }
  void initialize ( VelocityDiscreteFunctionType &velocity, PressureDiscreteFunctionType &pressure )
  {

    Dune::Fem::LagrangeInterpolation<InitialVelocityFunctionType, VelocityDiscreteFunctionType > mainVelocityInterpolation;
    mainVelocityInterpolation(mainModel_.initialVelocityFunction(), velocity );
      Dune::Fem::LagrangeInterpolation<InitialPressureFunctionType, PressureDiscreteFunctionType > mainPressureInterpolation;
     mainPressureInterpolation( mainModel_.initialPressureFunction(), pressure );
  }

  //! setup the right hand side
  void prepare( VelocityDiscreteFunctionType &velocity )
  {
    // set boundary values for velocity
    mainOperator_.prepare( mainModel_.dirichletBoundary(), velocity );

    // assemble rhs
    assembleRHS ( mainModel_, mainModel_.rightHandSide(), mainModel_.neumanBoundary(), rhsU_ );

  }

  void solve ( VelocityDiscreteFunctionType &velocity, PressureDiscreteFunctionType &pressure, bool assemble )
  {
    //! [Solve the system]
    if( assemble )
    {
      // assemble linear operator (i.e. setup matrix)
      mainOperator_.jacobian( velocity , mainLinearOperator_ );
      explicitMainOperator_.jacobian( velocity , explicitMainLinearOperator_ );
      gradOperator_.jacobian( pressure , gradLinearOperator_ );
      divOperator_.jacobian( velocity  , divLinearOperator_ );
      massOperator_.jacobian( pressure , massLinearOperator_ );
      precondOperator_.jacobian( pressure , precondLinearOperator_ );
    }

    MainLinearInverseOperatorType invMainOp( mainLinearOperator_, solverEps_, solverEps_ );
    MassLinearInverseOperatorType invMassOp( massLinearOperator_, solverEps_, solverEps_ );
    PrecondLinearInverseOperatorType invPrecondOp( precondLinearOperator_, solverEps_, solverEps_ );

    //Prepare the RHS f^* + Au^n - bu^n
    explicitMainOperator_(velocity, xi_);
    rhsU_+=xi_; xi_.clear();
    transportOperator_(velocity, xi_);
    rhsU_-=xi_; xi_.clear();

    //adjust for Schur complement
    gradOperator_( pressure, xi_ );
    rhsU_ -= xi_;
    //apply BCs
    mainOperator_.prepare( velocity, rhsU_ );

    //
    //
    //
    // A   := mainOperator
    // B   := divOperator
    // B^T := gradOperator
    //
    // b_1 should be rhsU but from prepare step:
    // b_1 - B^Tx_2 :=rhsU
    //
    // b_2 (=0 here) := rhsP
    //
    //Method
    //Schur: (A B^T)(x_1)=(b_1           )
    //       (0 -S )(x_2) (b_2-BA^{-1}b_1) where S=BA^{-1}B^T
    //
    //
    // We wish to solve      Sx_2 = BA^{-1}b_1 - b_2, for x_2
    // and perform recovery  Ax_1 = b_1 - B^Tx_2,     for x_1
    //
    // FIRST calc residual of Schur complement
    // r := BA^{-1}b_1 - b_2 - Sx_2 = A^{-1}(b_1-B^Tx_2) -b_2
    //
    // Choose d_2 := r search direction

    //  gradOperator_(pressure, xi_);
    //  rhsU_ -= xi_;


    invMainOp( rhsU_, velocity );//            solve Au =(b_1-B^Tx_2) => u = A^{-1}(b_1-B^Tx_2)
    divLinearOperator_( velocity, rhsP_ );//   set rhsP = BA^{-1}(b_1 -B^Tx_2)
    invMassOp( rhsP_, r_ );//                   solve mass*r =rhsP  => r=rhsP
    if (usePrecond_ && nu_>0.)
    {
      precond_.clear();                        //Preconditioning...
      invPrecondOp( rhsP_, precond_ );
      r_ *= mu_;
      r_.axpy(nu_,precond_);
    }
    d_.assign(r_);                             //set search direction d.
    double delta = r_.scalarProductDofs(rhsP_);// check preconditioning r^Tr
    assert( delta >= 0 );
    if ( delta < solverEps_*10. )             //delta = ||r|| check against tolerance
      return;

    // ENTER LOOP
    // compute a_2 := Sd_2 = BA^{-1}B^Td_2
    //   store d_1 := A^{-1}B^Td_2
    //
    // set scale factor rho = r^T r / d_2^T a_2
    //
    // UPDATES
    // set p = p - rho*d_2
    // set u = u - rho*d_1
    // set rnew = r - rho*a_2
    //
    // SEARCH CORRECTION (from Gramm Schmidt)
    // set lambda = -rnew^T rnew / r^Tr
    // set d_2 = r - lambda*d_2


    for (int m=0;m<100;++m) //any "sufficient" number of steps
    {
      xi_.clear();
      gradLinearOperator_(d_, rhsU_);                           // set rhsU =B^Td_2
      mainOperator_.prepare( mainModel_.zeroVelocity(), rhsU_ );// prepare Au = B^Td_2
      invMainOp( rhsU_, xi_ );                                  // solve Ax = B^Td_2  => store x= A^{-1}B^Td_2 (we call x "d_1")
      divLinearOperator_( xi_, rhsP_ );                         // set rhsP = Bx =BA^{-1}B^Td_2   (we call rhsP "a_2")
      double rho = delta / d_.scalarProductDofs(rhsP_);         // scale factor rho = r^Tr/(d_2^T S d_2) = delta/(d_2 . a_2)
      pressure.axpy(rho,d_);                                   //update pressure p=p-rho*d_2
      velocity.axpy(-rho,xi_);                                 //update velocity unew=u-rho*d_1
      divLinearOperator_( velocity, rhsP_ );                   // set rhsP = Bunew = Bu - rho*a_2
      invMassOp( rhsP_, r_ );                                   // solve mass*r=rhsP => r= r - rho*a_2
      if (usePrecond_ && nu_>0.)
      {
        precond_.clear();                                       // preconditioning...
        invPrecondOp( rhsP_, precond_ );
        r_ *= mu_;
        r_.axpy(nu_,precond_);
      }
      double oldDelta = delta;                                  // old ||r||
      delta = r_.scalarProductDofs(rhsP_);
      //std::cout << delta << std::endl;                    // new ||r||
      if ( delta < solverEps_*10. ) break;
      double gamma = delta/oldDelta;                            // "-1*"scale factor (we called lambda) to update search direction
      d_ *= gamma;                                              //set d_1 = (r_2 - lambda*d_1 =) r_2 + gamma*d_1
      d_ += r_;
    }
  }

protected:
  GridPartType &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with
  double mu_,nu_;

  const MainModelType mainModel_;
  const MainModelType explicitMainModel_;
  const GradModelType gradModel_;
  const DivergenceModelType divModel_;
  const MassModelType massModel_;
  const PrecondModelType precondModel_;
  const TransportModelType transportModel_;

  const VelocitySpaceType& velocitySpace_;
  const PressureSpaceType& pressureSpace_;
  VelocityDiscreteFunctionType rhsU_,xi_;
  PressureDiscreteFunctionType rhsP_,d_,r_,precond_;

  MainOperatorType mainOperator_,explicitMainOperator_;
  GradOperatorType gradOperator_;
  DivergenceOperatorType divOperator_;
  MassOperatorType massOperator_;
  PrecondOperatorType precondOperator_;
  TransportOperatorType transportOperator_;
  MainLinearOperatorType mainLinearOperator_,explicitMainLinearOperator_;
  GradLinearOperatorType gradLinearOperator_;
  DivergenceLinearOperatorType divLinearOperator_;
  MassLinearOperatorType massLinearOperator_;
  PrecondLinearOperatorType precondLinearOperator_;


  const double solverEps_ ; // eps for linear solver
  const bool usePrecond_;
};

#endif // end #if STOKES_UZAWASCHEME_HH
