// functions for optical model calculations

/**************************************/
/*      omcalc.cpp                    */
/**************************************/
int     omCalc      (double, Pdata *, ZAnumber *, double, double *, CrossSection *);

/**************************************/
/*      cccalc.cpp                    */
/**************************************/
int     ccCalc      (double, double, Pdata *, ZAnumber *, double, Direct *, double **, CrossSection *);
int     ccPreset    (double, double, Pdata *, ZAnumber *, double, Direct *);
void    ccMain      (double, double, Particle, ZAnumber *, double **, CrossSection *);
void    ccCleanUp   ();
int     ccHermiteMatrix (const int, const int);


/**************************************/
/*      dwcalc.cpp                    */
/**************************************/
int     dwbaCalc    (double, double, Pdata *, ZAnumber *, double, Direct *, CrossSection *);

