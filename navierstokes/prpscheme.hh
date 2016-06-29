#ifndef NAVIERSTOKES_PRPSCHEME_HH
#define NAVIERSTOKES_PRPSCHEME_HH

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
//#include "temporalprobleminterface.hh"
#include "navierstokes.hh"

#include "burgersmodel.hh"

#include "rhs.hh"
#include "elliptic.hh"
#include "noslipconstraints.hh"

#include "space.hh"
template <class DiscreteFunction>
double hSemiInnerProduct(const DiscreteFunction &u,const DiscreteFunction &w)
{

  typedef typename DiscreteFunction :: DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  const DiscreteFunctionSpaceType &dfSpace(u.space());
  //values of interest
  double semiInnerProduct = 0.0;
  //Loop over the elements
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
    {
      const ElementType &entity = *it;
      const GeometryType &geometry = entity.geometry();
      const LocalFunctionType uLocal = u.localFunction( entity );
      const LocalFunctionType wLocal = w.localFunction( entity );
      // obtain quadrature order
      const int quadOrder = uLocal.order();

      QuadratureType quadrature( entity, quadOrder );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

      JacobianRangeType du,dw;
      uLocal.jacobian( quadrature[ pt ], du );//2x2 matrix du_i/dx_j
      wLocal.jacobian( quadrature[ pt ], dw );//2x2 matrix dw_i/dx_j
      for (unsigned int i=0;i<GridPartType::dimensionworld;++i)
        {
          for (unsigned int j=0;j<GridPartType::dimensionworld;++j)
        {
          semiInnerProduct+=du[i][j]*dw[i][j]*weight; // A : B = sum_ij AijBij
        }
        }
    }
    }

  return u.gridPart().grid().comm().sum ( semiInnerProduct );
}

//------------------------------------------------------------//
//------------------------------------------------------------//
//------------------------------------------------------------//


template <class DiscreteFunction>
double energyFunctional(const DiscreteFunction &u,const double &alphaOne,const double &alphaTwo,const double &deltaT)
{

  typedef typename DiscreteFunction :: DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
  //typedef typename DiscreteFunctionType :: DomainType DomainType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
  //typedef typename GeometryType :: LocalCoordinate LocalCoordinateType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  //typedef Dune::Fem::TimeProviderBase TimeProviderType;

  const DiscreteFunctionSpaceType &dfSpace(u.space());
  //output vector
  //const double endTime  = Dune::Fem::Parameter::getValue< double >( "heat.endtime", 2.0 );
  //values of interest
  double energy = 0.0;
  double LTwoPart = 0.0;
  double HOnePart = 0.0;
  //Loop over the elements
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
    {
      const ElementType &entity = *it;
      const GeometryType &geometry = entity.geometry();
      const LocalFunctionType uLocal = u.localFunction( entity );
      // obtain quadrature order
      const int quadOrder = uLocal.order();

      QuadratureType quadrature( entity, quadOrder );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      LTwoPart=0.0;HOnePart=0.0;
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

      JacobianRangeType du,dw;
      uLocal.jacobian( quadrature[ pt ], du );//2x2 matrix du_i/dx_j
      for (unsigned int i=0;i<GridPartType::dimensionworld;++i)
        {
          for (unsigned int j=0;j<GridPartType::dimensionworld;++j)
        {
          HOnePart+=du[i][j]*du[i][j]; // A : B = sum_ij AijBij
        }
        }

      energy+=HOnePart*alphaTwo*0.5*weight;

      RangeType ru,rw;
      uLocal.evaluate( quadrature[ pt ], ru );//2-vector u_i
      for (unsigned int i=0;i<GridPartType::dimensionworld;++i)
        {
          LTwoPart+=ru[i]*ru[i]; // A . B = sum_i AiBi
        }
      energy+=LTwoPart*alphaOne*0.5*(1./deltaT)*weight;
    }
    }

  return u.gridPart().grid().comm().sum ( energy );
}


//------------------------------------------------------------//
//------------------------------------------------------------//
//------------------------------------------------------------//

template < class VelocitySpace, class PressureSpace >
class PRPScheme
{
public:
  typedef typename VelocitySpace::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;

  //Models we need to create
  typedef BurgersStateModel<GridPartType>      BurgersStateModelType;//This will construct the model L in disc notes
  typedef BurgersTransportModel<GridPartType>  BurgersTransportModelType;//This will construct the model L_B - L in disc notes
  typedef BurgersDescentModel<GridPartType>    BurgersDescentModelType;//This will construct the model D_J - L in disc notes
  typedef BurgersGradModel<GridPartType>       BurgersGradModelType;

  //Our function space
  typedef typename BurgersStateModelType::VelocityFunctionSpaceType VelocityFunctionSpaceType;
  typedef typename BurgersStateModelType::PressureFunctionSpaceType PressureFunctionSpaceType;
  //Our discrete function space
  typedef VelocitySpace VelocitySpaceType;
  typedef PressureSpace PressureSpaceType;
  //Our discrete functions
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< VelocitySpaceType > VelocityDiscreteFunctionType;
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< PressureSpaceType > PressureDiscreteFunctionType;
  typedef VelocityDiscreteFunctionType DiscreteFunctionType;
  //Invertable Operator type
  typedef Dune::Fem::ISTLLinearOperator< VelocityDiscreteFunctionType, VelocityDiscreteFunctionType > BurgersStateLinearOperatorType;
  typedef Dune::Fem::ISTLLinearOperator< PressureDiscreteFunctionType, VelocityDiscreteFunctionType > BurgersGradLinearOperatorType;
  //Inverse operators constructed from above
  typedef Dune::Fem::ISTLCGOp< VelocityDiscreteFunctionType, BurgersStateLinearOperatorType > BurgersStateLinearInverseOperatorType;

  typedef typename Dune::Fem::FunctionSpace<double,double,GridPartType::dimensionworld,GridPartType::dimensionworld+1> FullFunctionSpaceType;
  typedef NavierStokesProblemInterface<FullFunctionSpaceType> ProblemType ;
  //Define the linearisable/non-linearisable elliptic operator type + dirichlet constraints
  //typedef Dune::DirichletConstraints<BurgersStateModelType, VelocitySpaceType> StateConstraintsType;
  typedef Dune::DirichletConstraints<BurgersStateModelType, VelocitySpaceType> StateConstraintsType;
  typedef DifferentiableEllipticOperator< BurgersStateLinearOperatorType, BurgersStateModelType, StateConstraintsType > BurgersStateOperatorType;
  typedef EllipticOperator<  VelocityDiscreteFunctionType, VelocityDiscreteFunctionType, BurgersTransportModelType, NoConstraints > BurgersTransportOperatorType;
  typedef EllipticOperator<  VelocityDiscreteFunctionType, VelocityDiscreteFunctionType,  BurgersDescentModelType, NoConstraints > BurgersDescentOperatorType;
  typedef DifferentiableEllipticOperator< BurgersGradLinearOperatorType, BurgersGradModelType,NoConstraints > BurgersGradOperatorType;

  PRPScheme( const VelocitySpace &velocitySpace, const PressureSpace &pressureSpace,
             const ProblemType &problem,const double& alphaOne,const double& alphaTwo)
    : gridPart_( velocitySpace.gridPart() ),
      alphaOne_(alphaOne),
      alphaTwo_(alphaTwo),
      deltaT_( problem.deltaT() ),
      stateModel_( problem, gridPart_,alphaOne,alphaTwo,true),
      explicitStateModel_( problem, gridPart_,alphaOne,alphaTwo,false),
      gradModel_(gridPart_),
      transportModel_( problem,gridPart_ ),
      velocitySpace_( velocitySpace ),
      pressureSpace_( pressureSpace ),
      rhsU_( "rhsU", velocitySpace_ ),
      dummyOne_("dummyOne",velocitySpace_),
      dummyTwo_("dummyTwo",velocitySpace_),
      xi_("xi",velocitySpace_),
      descentModel_( problem,gridPart_,xi_),//requires memory location of xi_;
      g_("g",velocitySpace_),
      gdiff_("gdiff", velocitySpace_),
      d_("d",velocitySpace_),
      stateOperator_(stateModel_,velocitySpace_),
      explicitStateOperator_(explicitStateModel_,velocitySpace_),
      gradOperator_(gradModel_,velocitySpace_),
      transportOperator_( transportModel_, velocitySpace_ ),
      descentOperator_( descentModel_, velocitySpace_ ),
      // create linear operator (domainSpace,rangeSpace)
      stateLinearOperator_("assembled state operator",velocitySpace_,velocitySpace_),
      explicitStateLinearOperator_("assembled state operator",velocitySpace_,velocitySpace_),
      gradLinearOperator_("assembled gradient operator",pressureSpace_,velocitySpace_),
      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "stokes.solvereps", 1e-6 )),
      lineSearchAccept_(Dune::Fem::Parameter::getValue<double>("burgers.linesearchaccept", 0.001)),
      theta_(Dune::Fem::Parameter::getValue< double >("navierstokes.implicitfactor",0.585786))
  {
    // set all DoF to zero
    //solution_.clear();pressure_.clear();
    const std::string lineSearchMethods[] = { "backtracking","quadratic","cubic" };
    enum { backtracking, quadratic,cubic};
    const int lineSearchMethod = Dune::Fem::Parameter::getEnum("burgers.linesearchmethod", lineSearchMethods, backtracking);
    lineSearchMethod_=lineSearchMethod;
  }
  //! setup the right hand side
  void prepare( VelocityDiscreteFunctionType &solution )
  {
    // set boundary values for velocity
    stateOperator_.prepare( stateModel_.dirichletBoundary(), solution );
    // assemble rhs
    assembleRHS ( stateModel_, stateModel_.rightHandSide(), stateModel_.neumanBoundary(), rhsU_ );

    //NEED TO OBTAIN THE PRESSURE FROM STOKES.
  }

  void solve ( VelocityDiscreteFunctionType &solution, PressureDiscreteFunctionType &pressure, bool assemble )
  {

    VelocityDiscreteFunctionType solutiontemp(solution);//We wish to differentiate between solution at time n, and solutions in prior Conjugate gradient steps.

    //! [Solve the system]
    if( assemble )
    {
      // assemble linear operator (i.e. setup matrix)
      stateOperator_.jacobian( solutiontemp , stateLinearOperator_ );
    }

    BurgersStateLinearInverseOperatorType invStateOp( stateLinearOperator_, solverEps_, solverEps_ );

    //assemble RHS
    gradOperator_(pressure,dummyOne_);
    rhsU_-=dummyOne_; dummyOne_.clear();

    //-------------Polack-Ribiere-Polyak Scheme ---------------
    //u0 initial velocity in H^1 with bdry condition g
    //xi state equation solution in H^1_0
    //g direction of steepest slope, in H^1_0
    //phi finite element test function
    //alpha1,alpha2 constants >0;
    //First solve state equation:
    // alpha1*(xi,phi) + alpha2* ((xi,phi)) = alpha1*(u0,phi) +alpha2*((u0,phi))+((u0 cdot \nabla u0),phi)-(f,phi)
    //
    // L[u]  = state Operator
    //      := (alpha1*Mass + theta*alpha2*Stiffness)u
    //
    // LE[u] = explicit state Operator
    //      := (-alpha1*Mass + (1-theta)*alpha2*stiffness)u
    //
    // B[u]  = Burgers main Operator
    //      := (u\cdot \nabla u)*Mass
    //
    // D[u]  = Descent Operator
    //      := (xi(u)\otimes u)\colon Stiffness + \xi(u) \cdot \nabla u)* Mass
    //----------------------------------------------------------


    //----------------------Solve L[xi] = L[u0]+B[u0]-L[u^n]-RhsU
    stateLinearOperator_(solutiontemp,dummyOne_);                  // construct dummyOne = Lu
    transportOperator_(solutiontemp,dummyTwo_);                    // construct dummyTwo = Bu
    dummyTwo_ += dummyOne_; dummyTwo_ -= rhsU_;                               // set dummyTwo = Lu+Bu - f

    dummyOne_.clear();
    explicitStateOperator_(solution,dummyOne_);             // set dummyTwo = L[u0]+B[u0] + L'[u^n] - f
    dummyTwo_+= dummyOne_;///////

    stateOperator_.prepare(stateModel_.zeroVelocity(),dummyTwo_);
    invStateOp(dummyTwo_,xi_);                                    // solve Lxi =Lu+ Bu-f for xi (=> xi=L^{-1}(Lu+ Bu-f) )

    //----------------------Check condition:
    double delta = 0.0;
    delta = energyFunctional(xi_,alphaOne_,alphaTwo_,deltaT_);
    std::cout<<delta<<std::endl;
    if(delta < solverEps_/100.) return;

    //---------------------Solve for the steepest ascent direction:
    dummyOne_.clear();dummyTwo_.clear();
    stateLinearOperator_(xi_,dummyOne_);                           // dummyOne = Lxi
    descentOperator_(solutiontemp,dummyTwo_);                   // Construct dummyTwo = D[u]
    dummyTwo_ += dummyOne_;                                        // set dummyTwo = L[xi]+D[u]
    stateOperator_.prepare(stateModel_.zeroVelocity(),dummyTwo_);  // Solving problem with zero BCs
    invStateOp(dummyTwo_,g_);                                      // solve Lg= Lxi + Du for g


    //---------------------Define first search direction d = g
    d_.assign(g_);
    d_ *=-1.;     //steepest descent direction

    //Begin Loop
    double prpBetaDenominator;
    double prpBetaNumerator;
    double prpBeta;
    for (unsigned int m = 0; m < 100; ++m)
      {
          //--------------------Line search
    {
      switch(lineSearchMethod_)
        {
        case 0: //Backtracking
          {
        //Parameters
        double stepSize = 2.0;//step size double the actual start size
        double tolInit = lineSearchAccept_;//in (0,0.5)
        double tol;
        double reduction = 0.5; // how much to reduce step size if failed > 0.0
        double eValNew = 0.0;
        double eValOld=delta;
        VelocityDiscreteFunctionType solutiontmp(solutiontemp);
        VelocityDiscreteFunctionType xitmp(solutiontemp);

        tolInit*=g_.scalarProductDofs(d_);
        for(unsigned int m=0;m<20;++m)
          {
            stepSize=stepSize*std::pow(reduction,m+1);  // reduce stepsize  (1st loop reduction =1.0)
            tol=tolInit*stepSize;                       //=tol*g*d*stepSize cancels out first loop.

            solutiontmp.clear();
            solutiontmp.assign(solutiontemp);
            solutiontmp.axpy(stepSize,d_);              //=u+stepSize d

            //-------------------Solve state equation for proposed step
            xitmp.clear();dummyOne_.clear();dummyTwo_.clear();
            stateLinearOperator_(solutiontmp,dummyOne_);
            transportOperator_(solutiontmp,dummyTwo_);
            dummyTwo_ += dummyOne_; dummyTwo_-= rhsU_;

            dummyOne_.clear();
            explicitStateOperator_(solution,dummyOne_);// set dummyTwo = L[u0]+B[u0] - L[u^n] - f
            dummyTwo_+= dummyOne_;/////////

            stateOperator_.prepare(stateModel_.zeroVelocity(),dummyTwo_);
            invStateOp(dummyTwo_,xitmp);
            //calc new energy
            eValNew = energyFunctional(xitmp,alphaOne_,alphaTwo_,deltaT_);

            if(eValNew<=eValOld+tol)
              {
            xi_.assign(xitmp);
            solutiontemp.assign(solutiontmp);//update u^{n+1} =u^{n} + \lambda d^n
            break;
              }
          }
        break;
          }
        case 1: //quadratic line search
          {
        //Implement quadratic interpolation linesearch here
        break;
          }
        case 2: //cubic
          {
        //Implement cubic interpolation linesearch here
        break;
          }
        }
    }


    //-------------------Check condition: J(u^{n+1})<Tolerance=solverEps_*10
    delta = energyFunctional(xi_,alphaOne_,alphaTwo_,deltaT_);
    std::cout << delta << std::endl;
    if(delta < solverEps_/100.) break;

    //------------------Polak-Ribiere-Polyak beta & ascent direction
    gdiff_.assign(g_);//gdiff = g^n
    gdiff_ *= -1.0;

    prpBetaDenominator = 0.0;
    prpBetaDenominator += (1./deltaT_)*alphaOne_ * g_.scalarProductDofs(g_);// + int(g_i g_i )dx
    prpBetaDenominator += alphaTwo_ * hSemiInnerProduct(g_,g_);// + int (\nabla g_i \cdot \nabla g_i) dx

    dummyOne_.clear();dummyTwo_.clear();
    stateLinearOperator_(xi_,dummyOne_);                          // dummyOne = Lxi^{n+1}
    descentOperator_(solutiontemp,dummyTwo_);                  // Construct dummyTwo =Du^{n+1}
    dummyTwo_ += dummyOne_;
    stateOperator_.prepare(stateModel_.zeroVelocity(),dummyTwo_); // Solving problem with zero BCs
    invStateOp(dummyTwo_,g_);                                     // solve Lg^{n+1}= Du^{n+1} for g^{n+1}

    gdiff_+=g_;                                                               //gdiff = g^{n+1}-g^n
    prpBetaNumerator = 0.0;
    prpBetaNumerator += (1./deltaT_)*alphaOne_ * g_.scalarProductDofs(gdiff_);// + int(g_i gdiff_i )dx
    prpBetaNumerator += alphaTwo_ * hSemiInnerProduct(g_,gdiff_);             // + sum_iint (\nabla g_i \cdot \nabla gdiff_i) dx
    prpBeta=std::max(prpBetaNumerator/prpBetaDenominator,0.0);

    //--------------------update search direction
    d_ *= prpBeta;
         d_ -= g_; // d^{n+1} = -g^{n+1} + prpBeta*d^n       (recall g is gradient, -g descent dir)
      }
    solution.assign(solutiontemp);//update solution to scheme

  }
protected:
  GridPartType &gridPart_; // grid part(view), e.g. here the leaf grid the discrete space is build with
  double alphaOne_,alphaTwo_,deltaT_;

  const BurgersStateModelType stateModel_,explicitStateModel_;
  const BurgersGradModelType gradModel_;
  const BurgersTransportModelType transportModel_;


  const VelocitySpaceType& velocitySpace_;
  const PressureSpaceType& pressureSpace_;
  VelocityDiscreteFunctionType rhsU_,dummyOne_,dummyTwo_,xi_;
  const BurgersDescentModelType descentModel_;
  VelocityDiscreteFunctionType g_,gdiff_,d_;

  BurgersStateOperatorType stateOperator_,explicitStateOperator_;
  BurgersGradOperatorType gradOperator_;
  BurgersTransportOperatorType transportOperator_;
  BurgersDescentOperatorType descentOperator_;


  BurgersStateLinearOperatorType stateLinearOperator_,explicitStateLinearOperator_;
  BurgersGradLinearOperatorType gradLinearOperator_;
  const double solverEps_,lineSearchAccept_,theta_ ; // eps for linear solver
  int lineSearchMethod_;
};






#endif // end #if NAVIERSTOKES_PRPSCHEME_HH
