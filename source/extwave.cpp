/******************************************************************************/
/*  extwave.cpp                                                               */
/*        wave functions in free space                                        */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "optical.h"
#include "etc.h"
#include "terminate.h"

static inline double omGfunction(int, double, double, double, double);
static inline double omFfunction(int, double, double, double, double);
static inline double omDfunction(int, double, double, double, double);

/**********************************************************/
/*      I.A.Stegun and M.Abramowitz                       */
/*      Phys. Rev. 98, 1851(1955)                         */
/**********************************************************/
int omExternalFunction(int lmax, double rho_match, double coulomb, double coulomb_scat0,
                       Wavefunc *wfn)
{
  double g0,g1,p1,p2,p3;
  int    l;

  if(coulomb == 0.0){
    g0 = p1 = cos(rho_match);
    g1 = p2 = g0/rho_match+sin(rho_match);
    wfn->external[0] = complex<double>(g0,0.0);
    wfn->external[1] = complex<double>(g1,0.0);
  }
  else{
    omAsymptotic(rho_match, coulomb, coulomb_scat0, &p1, &p2);
    if((p1 == 0.0) && (p2 == 0.0)){
      message << "cannot calculate G-function for charged particle case";
//    cohNotice("omExternalFunction"); 
      wfn->external[0] = complex<double>(0.0,0.0);
      return(-1);
    }
    g0 = p1;
    g1 = p2 = ((coulomb+1.0/rho_match)*p1-p2)/(sqrt(1.0+coulomb*coulomb));
    wfn->external[0] = complex<double>(g0,0.0);
    wfn->external[1] = complex<double>(g1,0.0);
  }

  double a = (coulomb==0.0) ? 1.0 : cfmax(g0*g0,g1*g1);

  for(l=1 ; l<MAX_L-1 ; l++){
    p3 = omGfunction(l,coulomb,rho_match,p1,p2);
    wfn->external[l+1] = complex<double>(p3, 0.0);
    if( (a/(p3*p3)) < CRIT_LMAX && (l >= lmax) ) break;
    p1 = p2;  p2 = p3;
  }
  if(l == MAX_L-1){
    message << "maximum angular momentum " << l << " too large";
    cohNotice("omExternalFunction");
  }

  lmax = (l==MAX_L-1) ? l-1 : l;

  int  l0   = lmax+10;
  bool flag = false;
  do{
    flag = false;
    p1 = 0.0; p2 = 1.0e-36;
    for(l=l0 ; l>0 ; l--){
      p3 = omFfunction(l,coulomb,rho_match,p1,p2);
      if(l <= lmax+2) wfn->external[l-1].imag(p3);
      p1 = p2;  p2 = p3;
    }

    double f0 = wfn->external[0].imag();
    double f1 = wfn->external[1].imag();
    a =1.0/((f0*g1-g0*f1)*sqrt(1.0+coulomb*coulomb));

    for(l=0 ; l<=lmax+1 ; l++) wfn->external[l].imag(wfn->external[l].imag()*a);
    for(l=0 ; l<=lmax   ; l++){
      double wr = omDfunction(l,coulomb,rho_match,wfn->external[l  ].real(),
                                                  wfn->external[l+1].real() );
      double wi = omDfunction(l,coulomb,rho_match,wfn->external[l  ].imag(),
                                                  wfn->external[l+1].imag() );
      wfn->extderiv[l] = complex<double>(wr,wi);
      a = wfn->extderiv[l].imag() * wfn->external[l].real()
         -wfn->extderiv[l].real() * wfn->external[l].imag();

      if( fabs(a-1.0) > CRIT_WRONSK ){
        flag = true;
        l0 += 10;
        break;
      }
    }
    if(l0 > MAX_L0){
      message << "angular momentum " << l0 << " too large";
      cohTerminateCode("omExternalFunction");
    }
  }while(flag);

  return(lmax);
}


/**********************************************************/
/*      Coulomb Function -- G at Matching Radius          */
/**********************************************************/
double omGfunction(int l, double eta, double rho, double g1, double g2)
{
  double j  = (double)l;
  double j1 = (double)(l+1);
  return( ((j+j1)*(eta+j*j1/rho)*g2-j1*sqrt(j*j+eta*eta)*g1)
          /(j*sqrt(j1*j1+eta*eta)) );
}


/**********************************************************/
/*      Coulomb Function -- F at Matching Radius          */
/**********************************************************/
double omFfunction(int l, double eta, double rho, double g1, double g2)
{
  double j  = (double)l;
  double j1 = (double)(l+1);
  return( ((j+j1)*(eta+j*j1/rho)*g2-j*sqrt(j1*j1+eta*eta)*g1)
          /(j1*sqrt(j*j+eta*eta)) );
}


/**********************************************************/
/*      Derivarive of Coulomb Function at Matching Radius */
/**********************************************************/
double omDfunction(int l, double eta, double rho, double g1, double g2)
{
  double j1 = (double)(l+1);
  return( ((eta+j1*j1/rho)*g1-sqrt(j1*j1+eta*eta)*g2)/j1 );
}


/**********************************************************/
/*      Coulomb Function for Negative Energy              */
/**********************************************************/
int omExternalClosed(int lmax, double rho_match, double coulomb, Wavefunc *wfn)
{
  double p1, p2, p3;

  if(coulomb == 0.0){

    /*** recurrence form, f[n+1] = [2n+1]/x f[n] + f[n-1] */
    p1 = 1.0/rho_match;
    p2 = p1+1.0/(rho_match*rho_match);
    wfn->external[0] = complex<double>(p1,0.0);
    wfn->external[1] = complex<double>(p2,0.0);

    for(int l=1 ; l<lmax ; l++){
      p3 = (2.0*l+1.0)/rho_match*p2 + p1;
      wfn->external[l+1] = complex<double>(p3,0.0);
      p1 = p2;  p2 = p3;
    }

    /*** convert into k[n], by multiplying exp(-x) */
    double e = exp(-rho_match);
    for(int l=0 ; l<=lmax ; l++){
      wfn->external[l] *= e;
      wfn->external[l].imag(0.0);
    }

    /*** derivative, (xk[n])' = -n k[n] -x k[n-1] */
    wfn->extderiv[0].real(-exp(-rho_match));
    wfn->extderiv[0].imag(0.0);

    for(int l=1 ; l<=lmax ; l++){
      double wr = -l*wfn->external[l].real() - rho_match*wfn->external[l-1].real();
      wfn->extderiv[l] = complex<double>(wr,0.0);
    }
    /*** covert k[n] into x k[n] */
    for(int l=0 ; l<=lmax ; l++){
      double wr = rho_match*wfn->external[l].real();
      wfn->external[l].real(wr);
    }
  }
  else{
    /*** closed Coulomb function, F' at rho */
    double f1 = omAsymptoticClosed(rho_match, coulomb);

    int lp = lmax + 24 + (int)(5.0*coulomb);
    double x = 1.0;
    for(int l=lp ; l>=0 ; l--){
      double a = coulomb/(l+1.0);
      double b = a + (l+1.0)/rho_match;
      x  = -(a*a - 1.0)/(b + x) + b;
      if(l <= lmax) wfn->extderiv[l] = complex<double>(x, 0.0);
    }

    for(int l=0 ; l<=lmax ; l++){
      double g1 = wfn->extderiv[l].real();
      double d  = 1.0/sqrt(abs(g1 - f1));
      wfn->external[l] = complex<double>(d,d);
      wfn->extderiv[l] = complex<double>(g1*d,f1*d);

      double a = coulomb/(l+1.0);
      double b = a + (l+1.0)/rho_match;
      f1 = (a*a - 1.0)/(b - f1) - b;
    }
  }

  return(lmax);
}
