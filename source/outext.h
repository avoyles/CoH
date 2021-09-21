/*
   outext.h :
        prototype of functions for extra data output
 */

/**************************************/
/*      outext.cpp                    */
/**************************************/
void    extCumulativeLevels             (Nucleus *);
void    extContinuumDensity             (Nucleus *);
void    extFissionDensity               (const int);
void    extGammaBranch                  (const int);
void    extGammaLine                    (const double, double *, Nucleus *, const double);
void    extDecayWidth                   (const int, double *);
