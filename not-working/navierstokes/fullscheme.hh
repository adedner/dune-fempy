#ifndef STOKES_FULLSCHEME_HH
#define STOKES_FULLSCHEME_HH

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
#include "probleminterface.hh"

#include "stokesmodel.hh"

#include "space.hh"
#include "rhs.hh"
#include "elliptic.hh"
#include "noslipconstraints.hh"



template < class GridPart >
class FullScheme
{
public:
  typedef GridPart GridPartType;

  typedef typename Dune::Fem::FunctionSpace<double,double,GridPart::dimensionworld,GridPart::dimensionworld+1> FunctionSpaceType;
  typedef ProblemInterface< FunctionSpaceType > ProblemType ;
  typedef StokesModel<FunctionSpaceType, GridPart> ModelType;

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

#if !EQUAL_ORDER
  typedef Dune::Fem::FunctionSpace<double,double,GridType::dimensionworld,GridType::dimensionworld> VelocityFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace<double,double,GridType::dimensionworld,1> PressureFunctionSpaceType;
#if MINIELEMENT
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< VelocityFunctionSpaceType, GridPartType, POLORDER+1 > VelocitySpaceType;
#else
  typedef Dune::Fem::BubbleElementSpace< VelocityFunctionSpaceType, GridPartType > VelocitySpaceType;
#endif
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< PressureFunctionSpaceType, GridPartType, POLORDER > PressureSpaceType;
  typedef Dune::Fem::TupleDiscreteFunctionSpace< VelocitySpaceType, PressureSpaceType> DiscreteFunctionSpaceType;
#else
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;
#endif

  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::ISTLGMResOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;

  /*********************************************************/
  //! define no-slip constraints
  typedef Dune::NoSlipConstraints<ModelType, DiscreteFunctionSpaceType> ConstraintsType;

  //! define Stokes operator
  typedef DifferentiableEllipticOperator< LinearOperatorType, ModelType, ConstraintsType > EllipticOperatorType;

  FullScheme( GridPartType &gridPart,
              const ProblemType& problem,
              double mu, double nu)
    : model_( problem, gridPart, mu, nu ),
      gridPart_( gridPart ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      rhs_( "rhs", discreteSpace_ ),
      // the elliptic operator (implicit)
      implicitOperator_( model_, discreteSpace_ ),
      // create linear operator (domainSpace,rangeSpace)
      linearOperator_( "assembled elliptic operator", discreteSpace_, discreteSpace_ ),
      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "stokes.solvereps", 1e-8 ) )
  {
    // set all DoF to zero
    solution_.clear();
  }

  DiscreteFunctionType &solution()
  {
    return solution_;
  }
  const DiscreteFunctionType &solution() const
  {
    return solution_;
  }

  //! setup the right hand side
  void prepare()
  {
    // set boundary values for solution
    implicitOperator_.prepare( model_.dirichletBoundary(), solution_ );

    // assemble rhs
    assembleRHS ( model_, model_.rightHandSide(), model_.neumanBoundary(), rhs_ );

    // apply constraints, e.g. Dirichlet contraints, to the result
    implicitOperator_.prepare( solution_, rhs_ );
  }

  void solve ( bool assemble )
  {
    //! [Solve the system]
    if( assemble )
    {
      // assemble linear operator (i.e. setup matrix)
      implicitOperator_.jacobian( solution_ , linearOperator_ );
    }

    // inverse operator using linear operator
    LinearInverseOperatorType invOp( linearOperator_, solverEps_, solverEps_ );
    // solve system
    invOp( rhs_, solution_ );
    //! [Solve the system]
  }

protected:
  const ModelType model_;   // the mathematical model

  GridPartType  &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with

  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_;   // the unknown
  DiscreteFunctionType rhs_;        // the right hand side

  EllipticOperatorType implicitOperator_; // the implicit operator

  LinearOperatorType linearOperator_;  // the linear operator (i.e. jacobian of the implicit)

  const double solverEps_ ; // eps for linear solver
};

#endif // end #if STOKES_FULLSCHEME_HH
