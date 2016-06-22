#ifndef NAVIERSTOKES_PROBLEMS_HH
#define NAVIERSTOKES_PROBLEMS_HH

#include <cassert>
#include <cmath>

#include "temporalprobleminterface.hh"

//Framework class
template <class FunctionSpace>
class NavierStokesProblemInterface : public TemporalProblemInterface < FunctionSpace >
{
public:
  typedef FunctionSpace FunctionSpaceType;
  //typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
  typedef Dune::Fem::TimeProviderBase  TimeProviderType ;

  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };


public:

  NavierStokesProblemInterface( const Dune::Fem::TimeProviderBase &timeProvider)
    : BaseType(timeProvider){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
    phi=RangeType(0);
  }
    virtual void fnext(const DomainType& x,
                 RangeType& phi) const
  {
    phi=RangeType(0);
  }


  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    phi=RangeType(0);
  }
  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    ret=JacobianRangeType(0);
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
    return true;
  }

  virtual bool hasDirichletBoundary () const
  {
    return true;
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

protected:
  double viscosity_,sFactor_,bFactor_;
};

//------------------------------------------------------//
//------------------------------------------------------//
//------------------------------------------------------//

// Example flow
template <class FunctionSpace>
class ChannelFlow : public NavierStokesProblemInterface<FunctionSpace>
{
  typedef NavierStokesProblemInterface< FunctionSpace >  BaseType;
 typedef Dune::Fem::TimeProviderBase  TimeProviderType ;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };
public:

  ChannelFlow( const Dune::Fem::TimeProviderBase &timeProvider,double viscosity,double sFactor ,double bFactor)
    : BaseType(timeProvider),viscosity_(viscosity),sFactor_(sFactor),bFactor_(bFactor){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  double viscosity() const
  {
    return viscosity_;
  }
  double sFactor() const
  {
    return sFactor_;
  }
  double bFactor() const
  {
    return bFactor_;
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
       // velocity
    phi[0] = x[1] * ( 1.0 - x[ 1 ] );
    phi[1] = 0.;
    // pressure
    phi[ dimDomain ] = (-2.0 * x[0] + 2.0)*viscosity_;

  }
virtual void fnext(const DomainType& x,
                 RangeType& phi) const
  {
     phi = RangeType(0);

  }
  virtual void unext(const DomainType& x,
                 RangeType& phi) const
  {
    // velocity
    phi[0] = x[1] * ( 1.0 - x[ 1 ] );
    phi[1] = 0.;
    // pressure
    phi[ dimDomain ] = (-2.0 * x[0] + 2.0)*viscosity_;
  }
  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
     ret = JacobianRangeType(0);
  }
  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
                 RangeType& value) const
  {
    u(x,value);
  }
  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    return true;
  }

  virtual bool hasDirichletBoundary () const
  {
    return true;
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
  double viscosity_,sFactor_,bFactor_;
};

//---------------------------------------------------------//
//---------------------------------------------------------//
//---------------------------------------------------------//

template <class FunctionSpace>
class VorticityFlow : public NavierStokesProblemInterface<FunctionSpace>
{
  typedef NavierStokesProblemInterface< FunctionSpace >  BaseType;
 typedef Dune::Fem::TimeProviderBase  TimeProviderType ;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };
public:

  VorticityFlow( const Dune::Fem::TimeProviderBase &timeProvider,double viscosity,double sFactor ,double bFactor)
    : BaseType(timeProvider),viscosity_(viscosity),sFactor_(sFactor),bFactor_(bFactor){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  double viscosity() const
  {
    return viscosity_;
  }
  double sFactor() const
  {
    return sFactor_;
  }
  double bFactor() const
  {
    return bFactor_;
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
     phi = RangeType(0);
     phi *= viscosity_;
  }


  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    // velocity
    phi[0] = -cos(2.*M_PI*x[0])*sin(2.*M_PI*x[1])*exp(-8.*M_PI*M_PI*viscosity_*timeProvider().time());
    phi[1] = sin(2.*M_PI*x[0])*cos(2.*M_PI*x[1])*exp(-8.*M_PI*M_PI*viscosity_*timeProvider().time());
    // pressure
    phi[ dimDomain ] =-(1./4.)*(cos(4.*M_PI*x[0])+cos(4.*M_PI*x[1]))*exp(-16.*M_PI*M_PI*viscosity_*timeProvider().time());//1/visc = reynolds
  }
virtual void fnext(const DomainType& x,
                 RangeType& phi) const
  {
     phi = RangeType(0);
     phi *= viscosity_;

  }
  virtual void unext(const DomainType& x,
                 RangeType& phi) const
  {
     // velocity
    phi[0] = -cos(2.*M_PI*x[0])*sin(2.*M_PI*x[1])*exp(-8.*M_PI*M_PI*viscosity_*(timeProvider().time()+timeProvider().deltaT()));
    phi[1] = sin(2.*M_PI*x[0])*cos(2.*M_PI*x[1])*exp(-8.*M_PI*M_PI*viscosity_*(timeProvider().time()+timeProvider().deltaT()));
    // pressure
    phi[ dimDomain ] =-(1./4.)*(cos(4.*M_PI*x[0])+cos(4.*M_PI*x[1]))*exp(-16.*M_PI*M_PI*viscosity_*(timeProvider().time()+timeProvider().deltaT()));//1/visc = reynolds
  }
  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    ret=JacobianRangeType(0);
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
    return true;
  }

  virtual bool hasDirichletBoundary () const
  {
    return true;
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
  double viscosity_,sFactor_,bFactor_;
};

//---------------------------------------------------------//
//---------------------------------------------------------//
//---------------------------------------------------------//



template <class FunctionSpace>
class MovingPlate : public NavierStokesProblemInterface<FunctionSpace>
{
  typedef NavierStokesProblemInterface< FunctionSpace >  BaseType;
  typedef Dune::Fem::TimeProviderBase  TimeProviderType ;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };
public:

  MovingPlate( const Dune::Fem::TimeProviderBase &timeProvider,double viscosity,double sFactor ,double bFactor)
    : BaseType(timeProvider),viscosity_(viscosity),sFactor_(sFactor),bFactor_(bFactor){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  double viscosity() const
  {
    return viscosity_;
  }
  double sFactor() const
  {
    return sFactor_;
  }
  double bFactor() const
  {
    return bFactor_;
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
    phi[0]=0.0;
    double s = x[0]*x[0]*x[0]*x[0] -2.*x[0]*x[0]*x[0] + x[0]*x[0];//s
    double sx = 4.*x[0]*x[0]*x[0]-6.*x[0]*x[0] +2.*x[0];//s'
    double sxx = 12.*x[0]*x[0]-12.*x[0]+2.;//s''
    double sxxx = 24.*x[0]-12.;//s'''
    double S = (1./5.)*x[0]*x[0]*x[0]*x[0]*x[0] -(1./2.)*x[0]*x[0]*x[0]*x[0]+(1./3.)*x[0]*x[0]*x[0];//S=int(s)
    double Sone = s*sxx-sx*sx;//S1
    double Stwo = (1./2.)*s*s;
    double v = x[1]*x[1]*x[1]*x[1] -x[1]*x[1];
    double vy = 4.*x[1]*x[1]*x[1] - 2.*x[1];
    double vyy = 12.*x[1]*x[1] - 2.;
    double vyyy = 24.*x[1];
    double V = v*vyyy-vy*vyy;

    phi[1] =viscosity_*8.*(24.*S+2.*sx*vyy+sxxx*v)+64.*(Stwo*V-v*vy*Sone);

  }


  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {

    double s = x[0]*x[0]*x[0]*x[0] -2.*x[0]*x[0]*x[0] + x[0]*x[0];//s
    double sx = 4.*x[0]*x[0]*x[0]-6.*x[0]*x[0] +2.*x[0];//s'
    double S = (1./5.)*x[0]*x[0]*x[0]*x[0]*x[0] -(1./2.)*x[0]*x[0]*x[0]*x[0]+(1./3.)*x[0]*x[0]*x[0];//S=int(s)
    double Stwo = (1./2.)*s*s;
    double v = x[1]*x[1]*x[1]*x[1] -x[1]*x[1];
    double vy = 4.*x[1]*x[1]*x[1] - 2*x[1];
    double vyy = 12.*x[1]*x[1] - 2.;
    double vyyy = 24.*x[1];


    // If you want to smoothly increase the solution from 0 initial conditions.
    double smoother = 1.;
    // double smoother = std::min((timeProvider().time()/(10.*timeProvider().deltaT())),1.);
    phi[0]=smoother*8.*s*vy;
    phi[1]=smoother*(-8.*sx*v);
    phi[2]=smoother*(viscosity_*8.*(S*vyyy+sx*vy)+64.*(Stwo*(v*vyy-vy*vy)));
  }
  virtual void fnext(const DomainType& x,
             RangeType& phi) const
  {

      phi[0]=0.0;
    double s = x[0]*x[0]*x[0]*x[0] -2.*x[0]*x[0]*x[0] + x[0]*x[0];//s
    double sx = 4.*x[0]*x[0]*x[0]-6.*x[0]*x[0] +2.*x[0];//s'
    double sxx = 12.*x[0]*x[0]-12.*x[0]+2.;//s''
    double sxxx = 24.*x[0]-12.;//s'''
    double S = (1./5.)*x[0]*x[0]*x[0]*x[0]*x[0] -(1./2.)*x[0]*x[0]*x[0]*x[0]+(1./3.)*x[0]*x[0]*x[0];//S=int(s)
    double Sone = s*sxx-sx*sx;//S1
    double Stwo = (1./2.)*s*s;
    double v = x[1]*x[1]*x[1]*x[1] -x[1]*x[1];
    double vy = 4.*x[1]*x[1]*x[1] - 2.*x[1];
    double vyy = 12.*x[1]*x[1] - 2.;
    double vyyy = 24.*x[1];
    double V = v*vyyy-vy*vyy;

    phi[1] =viscosity_*8.*(24.*S+2.*sx*vyy+sxxx*v)+64.*(Stwo*V-v*vy*Sone);

  }
  virtual void unext(const DomainType& x,
             RangeType& phi) const
  {
    double s = x[0]*x[0]*x[0]*x[0] -2.*x[0]*x[0]*x[0] + x[0]*x[0];//s
    double sx = 4.*x[0]*x[0]*x[0]-6.*x[0]*x[0] +2.*x[0];//s'
    double S = (1./5.)*x[0]*x[0]*x[0]*x[0]*x[0] -(1./2.)*x[0]*x[0]*x[0]*x[0]+(1./3.)*x[0]*x[0]*x[0];//S=int(s)
    double Stwo = (1./2.)*s*s;
    double v = x[1]*x[1]*x[1]*x[1] -x[1]*x[1];
    double vy = 4.*x[1]*x[1]*x[1] - 2*x[1];
    double vyy = 12.*x[1]*x[1] - 2.;
    double vyyy = 24.*x[1];

    phi[0]=8.*s*vy;
    phi[1]=-8.*sx*v;
    phi[2]= viscosity_*8.*(S*vyyy+sx*vy)+64.*(Stwo*(v*vyy-vy*vy));
  }
  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    ret = JacobianRangeType(0);
  }
  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
                 RangeType& value) const
  {
    value=RangeType(0);
    if (x[1]>(1.-1e-10))
      {
    value[0] = 16.*(x[0]*x[0]*x[0]*x[0] - 2.*x[0]*x[0]*x[0] + x[0]*x[0]);
    value[1]=0.;
      }
  }
  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    return true;
  }
  virtual bool hasDirichletBoundary () const
  {
    return true;
  }
  bool hasNeumanBoundary () const
  {
    return false;
  }
  virtual void n(const DomainType& x,
                 RangeType& value) const
  {
    value=RangeType(0);
  }
  private:
  double viscosity_,sFactor_,bFactor_;
};

template <class FunctionSpace>
class CouetteFlow : public NavierStokesProblemInterface<FunctionSpace>
{
  typedef NavierStokesProblemInterface< FunctionSpace >  BaseType;
  typedef Dune::Fem::TimeProviderBase  TimeProviderType ;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };
public:

  CouetteFlow( const Dune::Fem::TimeProviderBase &timeProvider,double viscosity,double sFactor ,double bFactor)
    : BaseType(timeProvider),viscosity_(viscosity),sFactor_(sFactor),bFactor_(bFactor){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  double viscosity() const
  {
    return viscosity_;
  }
  double sFactor() const
  {
    return sFactor_;
  }
  double bFactor() const
  {
    return bFactor_;
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
    phi=RangeType(0);
    double v=0.5;
          phi[0] = v*x[1];
    phi[1]=0.;
  }
virtual void fnext(const DomainType& x,
                 RangeType& phi) const
  {
     phi = RangeType(0);
  }
  virtual void unext(const DomainType& x,
                 RangeType& phi) const
  {
    phi=RangeType(0);
  }
  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    ret = JacobianRangeType(0);
  }
  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
                 RangeType& value) const
  {
    value=RangeType(0);
    double v=0.5;
    if(((x[0]>0.)&&(x[0]<1))&&(x[1]>(1.-1e-10)))
      {
    value[0] = v;
    // value[1] = x[1] * (1.0-x[1]);
    value[1]=0.;
      }
    else if ((x[0]>(1.-1e-10))||(x[0]<1e-10))
      {
    value[0] = v*x[1];
    value[1]=0.;
      }
  }
  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    return true;
    /* if(((x[1]>0.)&&(x[1]<1))&&(x[0]<=1e-10)){return true;}
       else{return false;}*/
  }

  virtual bool hasDirichletBoundary () const
  {
    return true;
  }
  bool hasNeumanBoundary () const
  {
    return false;
  }
  virtual void n(const DomainType& x,
                 RangeType& value) const
  {
    value=RangeType(0);
  }
  private:
  double viscosity_,sFactor_,bFactor_;
};

//------------------------------------------------------//
//------------------------------------------------------//
//------------------------------------------------------//

// Karman vortex street flow
//DOMAIN [-1,5] x [-1,1]
template <class FunctionSpace>
class KarmanVortexStreet : public NavierStokesProblemInterface<FunctionSpace>
{
  typedef NavierStokesProblemInterface< FunctionSpace >  BaseType;
  typedef Dune::Fem::TimeProviderBase  TimeProviderType ;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };
public:

  KarmanVortexStreet( const Dune::Fem::TimeProviderBase &timeProvider,double viscosity,double sFactor ,double bFactor)
    : BaseType(timeProvider),viscosity_(viscosity),sFactor_(sFactor),bFactor_(bFactor){}

  const TimeProviderType& timeProvider() const
  {
    return BaseType::timeProvider();
  }
  double viscosity() const
  {
    return viscosity_;
  }
  double sFactor() const
  {
    return sFactor_;
  }
  double bFactor() const
  {
    return bFactor_;
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
    phi = RangeType(0);
  }
  virtual void fnext(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType(0);

  }
  virtual void unext(const DomainType& x,
             RangeType& phi) const
  {
    phi = RangeType(0);
  }
  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
     ret = JacobianRangeType(0);
  }
  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
                 RangeType& value) const
  {
    value=RangeType(0);
    //  if ((x[0]<(-1+1e-8))&&(std::abs(x[1])<0.1))
    if (x[0]<(-1+1e-8))
      {
    value[0]=std::min(1.0,(((x[1]+1.)*(1.-x[1])*timeProvider().time())/(10*timeProvider().deltaT())));
      }
  }
  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    if (std::abs(x[1])>(1-1e-8)){return true;}
    else if ((std::abs(x[0])<0.5)&&(std::abs(x[1])<0.5)){return true;}
    else if (x[0]<(-1+1e-8)){return true;}
    else {return false;}
  }

  virtual bool hasDirichletBoundary () const
  {
    return true;
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
  double viscosity_,sFactor_,bFactor_;
};


#endif // #ifndef NAVIERSTOKES_PROBLEMS_HH
