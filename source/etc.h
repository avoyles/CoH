// miscellaneous functions, complex variable manipulation carried over from old C code

#ifndef __COMPLEX_H__
#define __COMPLEX_H__
#include <complex>
#endif


static inline double cfmin(double x, double y){ return( (x<y) ? x : y); }
static inline double cfmax(double x, double y){ return( (x>y) ? x : y); }


/**************************************/
/*      etc.cpp                       */
/**************************************/
double  laguerre              (int,    double, double);
double  gam                   (double);
double  loggamma              (double);
double  legendre              (int, double);
double  legendre1             (int, double);
double  dlegendre             (int, double);
double  assocLegendrePol      (int,int,double);
double  SphericalHarmonics    (int,int,double,double);
double  bessi2                (int, double);
double  bessk2                (int, double);
double  expint                (double);
complex<double> expintE1      (complex<double>);

