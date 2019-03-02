#ifndef BURGERS_PROBLEMS_HH
#define BURGERS_PROBLEMS_HH

#include <cassert>
#include <cmath>

#include "temporalprobleminterface.hh"

// Example channel flow
template <class FunctionSpace>
class InverseExponential: public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
 typedef Dune::Fem::TimeProviderBase  TimeProviderType ;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  InverseExponential( const Dune::Fem::TimeProviderBase &timeProvider,const double & mu, const double &nu,const double& alphaOne,const double &alphaTwo,const double &burgersTheta)
    : BaseType(timeProvider),mu_(mu), nu_(nu) ,alphaOne_(alphaOne),alphaTwo_(alphaTwo),burgersTheta_(burgersTheta){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  double alphaOne() const
  {
    return alphaOne_;
  }
  double alphaTwo() const
  {
    return alphaTwo_;
  }
 double burgersTheta() const
  {
    return burgersTheta_;
  }
  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
     phi = RangeType(0);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
     phi[0] = 0.75 - 0.25*1./(1+exp( ( -4.*x[0] +4.*x[1]-timeProvider().time() ) * 80./32. ));
     phi[1] = 0.75 + 0.25*1./(1+exp( ( -4.*x[0] +4.*x[1]-timeProvider().time() ) * 80./32. ));
     //phi[0] = 0.75 - 0.25*1./( 1+exp( (-4.*x[0] +4.*x[1])* 80./32) );
     // phi[1] = 0.75 + 0.25*1./( 1+exp( (-4.*x[0] +4.*x[1])* 80./32) );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    // fix later
    ret = JacobianRangeType(0);
  }

  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
                 RangeType& value) const
  {
    u(x, value);
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    //return ( std::abs( x[ 0 ] - 1.0 ) > 1e-8 );
    return true;
  }

  virtual bool hasDirichletBoundary () const
  {
    return true;
    //return false;
  }
  bool hasNeumanBoundary () const
  {
    return false;
  }
  virtual void n(const DomainType& x,
                 RangeType& value) const
  {
    value = RangeType(0);
  }

  private:
  double mu_,nu_,alphaOne_,alphaTwo_,burgersTheta_;
};


// Example trig
template <class FunctionSpace>
class Trigonometric: public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
 typedef Dune::Fem::TimeProviderBase  TimeProviderType ;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  Trigonometric( const Dune::Fem::TimeProviderBase &timeProvider,double mu, double nu,double alphaOne,double alphaTwo)
    : BaseType(timeProvider),mu_(mu), nu_(nu) ,alphaOne_(alphaOne),alphaTwo_(alphaTwo){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  double alphaOne() const
  {
    return alphaOne_;
  }
  double alphaTwo() const
  {
    return alphaTwo_;
  }
  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
     phi = RangeType(0);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {

    phi[0] = sin(M_PI*x[0])+cos(M_PI*x[1]);
    phi[1] = x[0] + x[1];

  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    // fix later
    ret = JacobianRangeType(0);
  }

  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
                 RangeType& value) const
  {
    u(x, value);
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    //return ( std::abs( x[ 0 ] - 1.0 ) > 1e-8 );
    return true;
  }

  virtual bool hasDirichletBoundary () const
  {
    return true;
    //return false;
  }
  bool hasNeumanBoundary () const
  {
    return true;
  }
  virtual void n(const DomainType& x,
                 RangeType& value) const
  {
    value = RangeType(0);
  }

  private:
  double mu_,nu_,alphaOne_,alphaTwo_;
};


// Example Linear
template <class FunctionSpace>
class Linear: public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
 typedef Dune::Fem::TimeProviderBase  TimeProviderType ;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  Linear( const Dune::Fem::TimeProviderBase &timeProvider,double mu, double nu,double alphaOne,double alphaTwo)
    : BaseType(timeProvider),mu_(mu), nu_(nu) ,alphaOne_(alphaOne),alphaTwo_(alphaTwo){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  double alphaOne() const
  {
    return alphaOne_;
  }
  double alphaTwo() const
  {
    return alphaTwo_;
  }
  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
     phi = RangeType(0);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {

    phi[0] = (x[0]+x[1]-2.*x[0]*timeProvider().time())/(1-2.*timeProvider().time()*timeProvider().time());
    phi[1] = (x[0]-x[1]-2.*x[1]*timeProvider().time())/(1-2.*timeProvider().time()*timeProvider().time());

  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    // fix later
    ret = JacobianRangeType(0);
  }

  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
                 RangeType& value) const
  {
    u(x, value);
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    //return ( std::abs( x[ 0 ] - 1.0 ) > 1e-8 );
    return true;
  }

  virtual bool hasDirichletBoundary () const
  {
    return true;
    //return false;
  }
  bool hasNeumanBoundary () const
  {
    return true;
  }
  virtual void n(const DomainType& x,
                 RangeType& value) const
  {
    value = RangeType(0);
  }

  private:
  double mu_,nu_,alphaOne_,alphaTwo_;
};

// Example Linear
template <class FunctionSpace>
class Gaussian: public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
 typedef Dune::Fem::TimeProviderBase  TimeProviderType ;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  Gaussian( const Dune::Fem::TimeProviderBase &timeProvider,double mu, double nu,double alphaOne,double alphaTwo)
    : BaseType(timeProvider),mu_(mu), nu_(nu) ,alphaOne_(alphaOne),alphaTwo_(alphaTwo){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  double alphaOne() const
  {
    return alphaOne_;
  }
  double alphaTwo() const
  {
    return alphaTwo_;
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType(0);
    if (((x[0]<0.51)&&(x[1]<0.51))&&((x[0]>0.49)&&(x[1]>0.49)))
      {
        phi[0]=1.0;
        phi[1]=1.0;
      }
  }

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType(0);
    // phi[0] = 1./sqrt(2.*M_PI)*exp((-(x[0]+x[1]-0.5)*(x[0]+x[1]-0.5))/(4.*timeProvider().time()));
    // phi[1] = 0.0;
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    // fix later
    ret = JacobianRangeType(0);
  }

  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
                 RangeType& value) const
  {
    u(x, value);
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    //return ( std::abs( x[ 0 ] - 1.0 ) > 1e-8 );
    return true;
  }

  virtual bool hasDirichletBoundary () const
  {
    return true;
    //return false;
  }
  bool hasNeumanBoundary () const
  {
    return true;
  }
  virtual void n(const DomainType& x,
                 RangeType& value) const
  {
    value = RangeType(0);
  }

  private:
  double mu_,nu_,alphaOne_,alphaTwo_;
};

#endif // #ifndef NAVIERSTOKES_PROBLEMS_HH
