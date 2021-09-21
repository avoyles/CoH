/*
   kcksys.h : 
        define level density parameter file location
        prototype of function to read level density parameters
 */


#ifndef __DIR_H__
#define __DIR_H_
#include "dir.h"
#endif

/**************************************/
/*      kcksyst.cpp                   */
/**************************************/
int     kckDataRead                     (ZAnumber *, LevelDensity *);
double  kckAsymptoticLevelDensity       (double);
double  kckSpinCutoff                   (double);
double  kckTemperature                  (double, double);
double  kckE0                           (double, double, double);
