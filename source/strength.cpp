/******************************************************************************/
/*  strength.cpp                                                              */
/*        strength function calculation                                       */
/******************************************************************************/

#include <iostream>
#include <cmath>

using namespace std;

#include "optical.h"
#include "omoutput.h"
#include "etc.h"


static const int MAX_STRENGTH_FUNCTION = 4;   /* = s, p, d, and f waves */

static complex<double> stfRmatrix   (double, complex<double>, complex<double>, complex<double>, complex<double>);
static double  stfRadius    (int);

static double pex = 0.0;


/**********************************************************/
/*      Strength Function calculation                     */
/*         Theory of Neutron Resonance Cross Section for  */
/*         Safety Applications                            */
/*                        F.H. Froehner, KfK 5073 (1992)  */
/**********************************************************/
void omStrengthFunction(int lmax0, int anum, double k2, double ein,
                        Wavefunc *wfn, complex<double> *smt)
{
  double  x,sr,si;
  complex<double> b,s,stf[MAX_STRENGTH_FUNCTION],rmt[MAX_STRENGTH_FUNCTION];

  double rm = stfRadius(anum);
  double a  = rm*sqrt(k2);
 
  int lmax1 = omExternalFunction(0,a,0.0,0.0,wfn);
  int lmax  = min(MAX_STRENGTH_FUNCTION-1,min(lmax0,lmax1));

  b = complex<double>(0.0,0.0);
  rmt[0]    = stfRmatrix(a,smt[0],b,wfn->external[0],wfn->extderiv[0]);

  double  p = 2*pex*sqrt(1.0e-06/ein)/PI;

  sr = ((x = 1.0-rmt[0].real()) < 0.0) ? 0.0 : rm*x;
  si = p*rmt[0].imag();
  stf[0] = complex<double>(sr, si);

  for(int l=1 ; l<=lmax ; l++){
    /*** j-averaged S-matrix */
    s = complex<double>((l+1.0)/(2.0*l+1.0),0.0) *smt[3*l  ]
       +complex<double>((l    )/(2.0*l+1.0),0.0) *smt[3*l+1];
    /*** boundary condition taken to be -l: Wigner-Eisenbad */
    b = complex<double>(-l,0.0);

    /*** energy-averaged R-matrix element, <R> */
    rmt[l]  = stfRmatrix(a,s,b,wfn->external[l],wfn->extderiv[l]);

    /*** <R> = R^\infty + i \pi s */
    /*** R' = a[1 - (2l+1)R^\infty]^{1/(2l+1)} */
    sr = ((x=1.0-(2*l+1)*rmt[l].real()) < 0.0) ? 0.0 : rm*pow(x,1.0/(2*l+1.0));
    /*** pole strength times 2 P0 sqrt(1eV/En) / pi */
    si = p*rmt[l].imag();

    stf[l] = complex<double>(sr, si);
  }
  outRmatrix(lmax,stf,rmt);
}


/**********************************************************/
/*      R-matrix from S-matrix                            */
/**********************************************************/
complex<double> stfRmatrix(double a, complex<double> s, complex<double> b, complex<double> w, complex<double> dw)
{
  /*** phase = tan^{-1} F/G */
  double p =  atan(w.imag()/w.real());

  /*** D = a x (G'+iF')/(G+iF) = L */
  complex<double> d =  dw/w * a;

  /*** for s-wave, save imag(L) = P */
  if(b.real() == 0) pex = d.imag();

  /*** L0 = L - B */
  complex<double> d0 = d - b;
  complex<double> e  = complex<double>(cos(2*p),-sin(2*p));
  /*** exp(-2ip) - S */
  complex<double> e1 = e - s;
  /*** (L -B ) (exp(-2ip) - S) */
  complex<double> c  = d0*e1;
  /*** (L - B) (exp(-2ip) - S) - 2iP exp(-2ip) */
  complex<double> e2(c.real()+2*d.imag()*e.imag(), c.imag()-2*d.imag()*e.real());
  /*** R = exp(-2ip) - S / {(L - B) (exp(-2ip) - S) - 2iP exp(-2ip)} */
  complex<double> r  = e1/e2;

  return(r);
}


/**********************************************************/
/*      Define Radius from A                              */
/**********************************************************/
inline double stfRadius(int a)
{
   return( pow((double)a,1.0/3.0) * 1.23 + 0.80 );
}
