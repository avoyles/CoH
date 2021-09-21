#ifndef __OUTPREP_H__
#define __OUTPREP_H__
#include "outprep.h"
#endif

#ifndef __BETACHAIN_H__
#define __BETACHAIN_H__
#include "betachain.h"
#endif


/**************************************/
/*      beohoutput.cpp                */
/**************************************/
void    outBanner                  (void);
void    outZA                      (ZAnumber *);
void    outTitle                   (char *);
void    outSystem                  (CalcMode, System *);
void    outFissionFragment         (const int, FFragData *);
void    outTargetState             (NuclearStructure *, const int);
void    outCompound                (const int, Pdata *);
void    outLevelDensity            (const int, double);
void    outGDR                     (GDR *);
void    outFissionBarrier          (int);
void    outBeta                    (ZAnumber *, Beta *, Beta *);
void    outBetaProfile             (BetaProfile *);
void    outGSProduction            (const int);
void    outFission                 (const int);
void    outSpectrum                (const bool, const double, double **);
void    outSpectrumSum             (const bool, const double, double **);
void    outSpectrumLab             (const int, const double, double *);
void    outGammaCascade            (Nucleus *);
void    outTotalResidual           (CumulativeResidualProduct *);
void    outIsomericRatio           (CumulativeResidualProduct *);
void    outFissionProductYield     (const int, Isotope *);
void    outFissionProductChainYield(const int, Isotope *);
void    outDelayedNeutronYield     (const int, Isotope *);

