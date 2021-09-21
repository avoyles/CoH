/******************************************************************************/
/*  dsdoverlap.cpp                                                            */
/*        overlap integral <R|f|x>                                            */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "structur.h"
#include "optical.h"
#include "dsd.h"

static complex<double> dsdVibCouplingFactor (int, double, double, double);

/**********************************************************/
/*      Direct/Semidirect Formfactors                     */
/**********************************************************/
void    dsdFormfactor(int c, Optical *omp, Potential *pot,
                      double v1, double w1, complex<double> *f1, complex<double> *f2)
{
  for(int i=0 ; i<=pot->n_match ; i++){
    double r=(double)i*(pot->width);

    /*** direct formfactor,  <R|r|X> */
    f1[i] = complex<double>(r,0.0);

    /*** semidirect formfactor,  <R|h(r)|X> */
    complex<double> h = dsdVibCouplingFactor(c, r,omp->R0 ,omp->a0);
    f2[i] = complex<double>(v1*h.real(),w1*h.imag());
  }
}


/**********************************************************/
/*      Semidirect Coupling Factor                        */
/**********************************************************/
complex<double> dsdVibCouplingFactor(int c, double x, double r, double a)
{
  complex<double> f;

  double z = exp((x-r)/a);
  switch(c){
  case 0:    // J.P. Boisson and S.Jang,  df(r)/dr
             // Nucl. Phys. A189, 334 (1972)
    f = complex<double>(z/((1+z)*(1+z)),0.0);
    break;
  case 1:    // H. Kitazawa, rf(r)
             // Nucl. Phys. A307, 1 (1978)
    f = complex<double>(x/(1+z), 0.0);
    break;
  case 2:    // M. Potokar, rf(r)-i W 4a df(r)/dr
             // Phys. Lett., 46B, 346 (1973)
    f = complex<double>(x/(1+z), x*z/((1+z)*(1+z)) *4);
    break;
  default:
    f = complex<double>(0.0, 0.0);
  }

  return(f);
}


/**********************************************************/
/*      Overlap Integral                                  */
/*      Distorted Wave and Bound Wave                     */
/**********************************************************/
complex<double> dsdRadialIntegral(int nint, double width, complex<double> *form,
                          double *bw, complex<double> *wfn)
{
  complex<double> sum;

  double f0=width    /3.0;
  double f1=width*2.0/3.0;
  double f2=width*4.0/3.0;

  sum = wfn[nint] * form[nint] * bw[nint] * f0;

  for(int i=1 ; i<nint ; i++){
    sum += wfn[i] * form[i] * bw[i] * ((i%2==0) ? f1 : f2);
  }

  return(sum);
}
