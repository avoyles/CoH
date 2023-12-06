/******************************************************************************/
/*  fispenetration.cpp                                                        */
/*        fission penetration factor by directly solving Schroedinger Eq.     */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "structur.h"
#include "statmodel.h"
#include "levden.h"
#include "nucleus.h"
#include "fpotential.h"
#include "terminate.h"


static inline complex<double>
lagrange(double h, complex<double> a,complex<double> b, complex<double> c,
                   complex<double> e,complex<double> f, complex<double> g)
{
  return( ((g-a)/60.0+0.15*(b-f)+0.75*(e-c))/h );
}

static complex<double> cnrm(0.0,0.0);

/**********************************************************/
/*      One-Dimensional Potential Energy                  */
/**********************************************************/
int     fpotPotentialEnergy(const int nx, double mu, double dx,
                            complex<double> *vpot, Fission *fb)
{
  const double beta0 = 0.05;
  double xbeta[5], xjunc[4], c[5], vhigh[5], homeg[5];

  /*** if triple-humped, there will be five regions, otherwise three */
  int nregion = 5;
  if(fb->barrier[2].height == 0.0) nregion = 3;

  vhigh[0] = fb->barrier[0].height;
  homeg[0] = fb->barrier[0].curvature;

  vhigh[1] = fb->potwell[0].height;
  homeg[1] = fb->potwell[0].curvature;

  vhigh[2] = fb->barrier[1].height;
  homeg[2] = fb->barrier[1].curvature;

  if(nregion == 5){
    vhigh[3] = fb->potwell[1].height;
    homeg[3] = fb->potwell[1].curvature;

    vhigh[4] = fb->barrier[2].height;
    homeg[4] = fb->barrier[2].curvature;
  }
  else{
    vhigh[3] = vhigh[4] = homeg[3] = homeg[4] = 0.0;
  }
  
  for(int i=0 ; i<5 ; i++) c[i] = mu * homeg[i] * homeg[i];


  /*** xbeta is the center of each curves
       first point, set zero at beta = 0, and shift by beta0 */
  xbeta[0] = sqrt(2.0 * vhigh[0] / c[0]) + beta0;
  for(int i=1 ; i<nregion ; i++){
    double dv = abs(vhigh[i-1] - vhigh[i]);
    xbeta[i] = xbeta[i-1] + sqrt(2.0 * dv * (c[i-1] + c[i])/(c[i-1] * c[i]));
  }

  /*** calculate conjunction points */
  for(int i=0 ; i<4 ; i++) xjunc[i] = 0.0;
  for(int i=0 ; i<nregion-1 ; i++){
    xjunc[i] = (c[i]*xbeta[i] + c[i+1]*xbeta[i+1])/(c[i] + c[i+1]);
  }


  /*** check the connection point order */
  bool connect = false;
  if(   (xbeta[0] < xjunc[0]) && (xjunc[0] < xbeta[1])
     && (xbeta[1] < xjunc[1]) && (xjunc[1] < xbeta[2]) ) connect = true;

  if(connect && (nregion == 5)){
    if(   (xbeta[2] < xjunc[2]) && (xjunc[2] < xbeta[3])
       && (xbeta[3] < xjunc[3]) && (xjunc[3] < xbeta[4]) ) connect = true;
  }
  if(!connect){
    message << "paraboras cannot be connected : ";
    message << " " << xbeta[0] <<" "<< xjunc[0] << " " << xbeta[1] << " " << xjunc[1];
    message << " " << xbeta[2] <<" "<< xjunc[2] << " " << xbeta[3] << " " << xjunc[3];
    message << " " << xbeta[4];
    cohTerminateCode("fpotPotentialEnergy");
  }

  /*** initialize */
  for(int i=0 ; i<nx ; i++) vpot[i] = complex<double>(0.0,0.0);


  /*** shift potential, and set zero below i0 */
  int nm = 0, i0 = 0;
  for(int i=i0 ; i<nx ; i++){
    double vx = i * dx;
    double yr = vhigh[0] - 0.5*c[0]*(vx - xbeta[0])*(vx - xbeta[0]);
    double yi = 0.0;
    vpot[i] = complex<double>(yr,yi);
    if(vpot[i].real() < 0.0){
      vpot[i] = complex<double>(0.0,0.0);
    }else{
      i0 = i;
      break;
    }
  }


  /*** calculate potential energy */
  for(int i=i0 ; i<nx ; i++){
    double vx = i * dx;

    int j = 0;
    for(j=0 ; j<nregion-1 ; j++) if(vx < xjunc[j]) break;

    double yr = vhigh[j] + ((j%2==0) ? -1.0 : 1.0)*0.5*c[j]*(vx - xbeta[j])*(vx - xbeta[j]);
    double yi = 0.0;

    /*** class-II and class-III well, if imaginary part is given */
    if( (j == 1) || (j == 3) ){
      double w = (j == 1) ? fb->potwell[0].absorb : fb->potwell[1].absorb;
      yi = yr - vhigh[j] - w;
      if(yi > 0.0) yi = 0.0;
    }

    /*** complex potential until Re V < 0 */
    vpot[i] = complex<double>(yr,yi);
    if(vpot[i].real() < 0.0){
      vpot[i] = complex<double>(0.0,0.0);
      nm = i;
      break;
    }
  }
  if(nm == 0) return(nm);

  /*** add extra points for matching calculation */
  nm += 20;
  if(nm > nx-1) nm = nx-1;

  return(nm);
}


/**********************************************************/
/*      One-Dimensional Potential Energy                  */
/**********************************************************/
void    fpotInternalFunction(const int nm, double ex, double mu, double dx, double wave,
                             double *tc, complex<double> *vpot)
{
  double h2 = dx*dx/12.0;
  double c  = 2.0*mu;

  /*** wave function integration using Fox-Goodwin method */
  complex<double> p0, p1, p2, q0(0.0,0.0), q1(0.0,0.0), q2;
  p0 = complex<double>(cos(-dx*wave),-sin(-dx*wave));
  p1 = complex<double>(cos(0.0),sin(0.0));

  q0 = complex<double>(-ex*c,0.0);
  q1 = complex<double>(-ex*c,0.0);

  for(int i=0 ; i<nm ; i++){
    q2 = (vpot[i] - ex)*c;
    p2 = ((2.0 + 10.0*h2*q1)*p1 - (1.0 - h2*q0)*p0) / (1.0 - h2*q2);

    q0 = q1;  q1 = q2;
    p0 = p1;  p1 = p2;
  }

  /*** extra points outside the potential for wave function matching
       using 7-point Lagrange derivative */
  complex<double> px[7];
  for(int i=nm ; i<nm+7 ; i++){
    q2 = - ex * c;
    p2 = ((2.0 + 10.0*h2*q1)*p1 - (1.0 - h2*q0)*p0) / (1.0 - h2*q2);

    q0 = q1;  q1 = q2;
    p0 = p1;  p1 = p2;
    px[i-nm] = p2;
  }
  complex<double> f = lagrange(dx,px[0],px[1],px[2],px[4],px[5],px[6]) / px[3];

  /*** matching point, and free particle wave functions */
  double xm = (nm+3) * dx * wave;

  complex<double> u1(cos(xm), sin(xm));
  complex<double> u2(cos(xm),-sin(xm));
  complex<double> du1(-wave*sin(xm), wave*cos(xm));
  complex<double> du2(-wave*sin(xm),-wave*cos(xm));

  /*** S-matrix */
  complex<double> smat = (f*u2 - du2)/(f*u1 - du1);

  /*** normalization factor */
  cnrm = (u2 - smat * u1) / px[3];

  tc[0] = 1.0 - (smat.real()*smat.real() + smat.imag()*smat.imag());
  tc[1] = norm(cnrm);
}


/**********************************************************/
/*      Normalized Wave Function                          */
/**********************************************************/
void    fpotInternalNormalized(const int nm, double ex, double mu, double dx, double wave,
                               double *tc, complex<double> *vpot)
{
  double h2 = dx*dx/12.0;
  double c  = 2.0*mu;

  /*** wave function integration using Fox-Goodwin method */
  complex<double> p0, p1, p2, q0(0.0,0.0), q1(0.0,0.0), q2, px, v;
  p0 = complex<double>(cos(-dx*wave),-sin(-dx*wave));
  p1 = complex<double>(cos(0.0),sin(0.0));

  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  for(int i=0 ; i<nm+7 ; i++){

    v = (i >= nm) ? complex<double>(0.0,0.0) : vpot[i];

    q2 = (v - ex)*c;
    p2 = ((2.0 + 10.0*h2*q1)*p1 - (1.0 - h2*q0)*p0) / (1.0 - h2*q2);

    q0 = q1;  q1 = q2;
    p0 = p1;  p1 = p2;

    px = p2 * cnrm;
    cout << setw(12) << i*dx;
    cout << setw(12) << v.real() << setw(12) << v.imag();
    cout << setw(12) << px.real() << setw(12) << px.imag() << endl;
  }
}
