/*
   fpotential.h : 
        prototype of functions for fission potential
 */


#ifndef __COMPLEX_H__
#define __COMPLEX_H__
#include <complex>
#endif

/**************************************/
/*      fpotential.cpp                */
/**************************************/
int     fpotPotentialEnergy  (const int, double, double, complex<double> *, Fission *);
void    fpotInternalFunction (const int, double, double, double, double, double *, complex<double> *);
void    fpotInternalNormalized  (const int, double, double, double, double, double *, complex<double> *);
