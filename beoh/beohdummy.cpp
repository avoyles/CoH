/******************************************************************************/
/*  beohdummy.cpp                                                             */
/*        this is a dummy file, including functions not used in BeoH          */
/******************************************************************************/

#include <iostream>

using namespace std;

#include "structur.h"
#include "statmodel.h"

void statWidthFluctuationSet(int a, int b, int c, int d, double e)
{ cerr << a << b << c << d << e; }

void statWidthFluctuationReset(double e)
{ cerr << e; }

void statWidthFluctuation(double a, double b, double c)
{ cerr << a << b << c; }

void wfcheck(){}

double  statMoldauer(int a, int b, int c, int d, double e, double f)
{ return(a+b+c+d+e+f); }
