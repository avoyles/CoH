// functions for fission neutron spectrum model calculations

const double MAX_SPECTRUM_ENERGY = 30.0;

/**************************************/
/*      fnscalc.cpp                   */
/**************************************/
void    fnsFissionNeutronModel        (System *, Pdata *, double *, FNSpec *);
void    fnsMadlandNixModel            (int, int, double *, double **, FChance *);
void    fnsInverseCrossSection        (Pdata *, ZAnumber *, ZAnumber *);
double  fnsAverageSpectrumEnergy      (const int, double *, double *);
void    fnsNormalizeSpectrum          (const int, double *, double *);
void    fnsNormalizeExcitationDist    (const int, double *, double *);


/**************************************/
/*      fnsprefis.cpp                 */
/**************************************/
double  fnsPreFissionSpectrum         (int, double *, double *, double *, Nucleus *);
