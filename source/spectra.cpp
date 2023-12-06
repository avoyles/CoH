/******************************************************************************/
/*  spectra.cpp                                                               */
/*        main particle emission spectra calculation                          */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"
#include "output.h"
#include "eclipse.h"
#include "global.h"
#include "terminate.h"

static double  specRenormReaction    ();

#undef DEBUG_POPCHECK
#ifdef DEBUG_POPCHECK
static void    specPopCheck          (Nucleus *);
#endif


/**********************************************************/
/*      Decay of each compound nucleus                    */
/*      ----                                              */
/*      Although this subroutine name does not follow     */
/*      the naming standard here, this was taken from     */
/*      the historic GNASH code.                          */
/**********************************************************/
void    spectra
( System        *sys,              // system parameters
  Pdata         *pdt,              // data for all emitting partiles
  Transmission  *tin,              // transmission coeff. for incident channel
  Transmission  **tc,              // transmission coeff. for continuum
  Transmission  **td,              // transmission coeff. for discrete levels
  double        **tg,              // transmission coeff. for gamma-ray
  Direct        *dir,              // direct reaction data
  Spectra       *spc)              // particle emission spectra for each channel
{
  /*** If preeq and DSD exist, renormalize total CN */
  double renorm = specRenormReaction();

  /*** Loop over parent compound nucleus */
  /***
       c0: index of the parent CN
       k0: index of the CN excitation. k0 changes from the highest
           excitation energy down to lust above the highest discrete 
           level of the parent CN. We assume that discrete states do
           not emit particles.
  */

  for(int c0=0 ; c0<sys->max_compound ; c0++){

    /*** Store particle transmission data into array */
    int k0 = 0;
    if(c0 != 0) statStoreContinuumTransmission(c0,ncl[c0].excitation[k0],pdt,tc);

    /*** Clear spectrum array */
    spc->memclear("cn");

    /*** Loop over parent excitation */
    for(k0=0 ; k0<ncl[c0].ncont ; k0++){

#ifdef DEBUG_POPCHECK
      if(k0 == 0) specPopCheck(&ncl[c0]);
#endif

      spc->memclear("dp");

      /*** Store gamma-ray transmission data into array */
      statStoreContinuumGammaTransmission(k0,tg,&ncl[c0]);

      /*** Initial binary reaction, including compound elastic */
      if((c0 == 0) && (k0 == 0) && (sys->incident.pid != gammaray)){

        /*** Compound elastic calculation */
        specCompoundElastic(sys,renorm,tin,tc,td,tg,dir,spc);

        /*** Extract primary gamma-ray from level population */
        spc->npg = specStorePrimaryGamma(&ncl[c0],spc->pg,&gml);

        /*** Output discrete level population before gamma cascading */
        if(prn.levelexcite) outDiscreteLevelPopulation(sys->ex_total,&ncl[0]);

        /*** Outupt primary gamma-ray lines */
        if(prn.spectra) outPrimaryGammaSpectrum(spc->npg,spc->pg,&ncl[0]);

#ifdef DEBUG_POPCHECK
        specPopCheck(&ncl[c0]);
#endif

      /*** General CN decay */
      }else{

        /*** Transmission coefficients for discrete levels */
        statStoreDiscreteTransmission(c0,ncl[c0].excitation[k0],pdt,td);

        /*** All residual CN decay */
        specCompoundDecay(c0,k0,tc,td,tg,spc);
      }
    }

    /*** Particle decay from levels if possible */
    specLevelDecay(c0,pdt,&ncl[c0],tc,td,spc);

    /*** Gamma-ray cascading  */
    specGammaCascade(spc->cn[0],&ncl[c0]);

    /*** Add ground state population */
    crx.prod[c0].xsec += ncl[c0].lpop[0];

    /*** Add to total particle emission spectra */
    specCumulativeSpectra(spc->getCsize(),spc->getNsize(),spc->cn,&ncl[c0]);

    /*** Output calculated spectra and populations */
    if(prn.spectra)     outSpectrum(0,spc->cn,&ncl[c0]);
    if(prn.xsection)    outSpectrumSum(0,spc->cn,&ncl[c0]);
    if(prn.gcascade)    outGammaCascade(c0,&ncl[c0]);
    if(prn.xsection)    outIsomerProduction(&ncl[c0]);
    if(prn.population)  outPopulation(&ncl[c0]);

    /*** Store discrete gamma lines */
    if(pex.gammaline || opt.finegammaspectrum) specStoreDiscreteGamma(&ncl[c0],&gml);
  }

  /*** Output decay widths of initial compound nucleus */
  if(pex.decaywidth) statDecayWidth(sys,pdt,tc,td,tg,spc);
}


/**********************************************************/
/*      Data Generation for Exclusive Spectrum Calc.      */
/*      ----                                              */
/*      Same calculation procedure as spectra routine.    */
/*      The output data are processed by ECLIPSE          */
/**********************************************************/
void    spectraExclusive(System *sys, Pdata *pdt, Transmission *tin,
                         Transmission **tc, Transmission **td,
                         double **tg, Direct *dir, Spectra *spc)
{
  /*** discrete population increment other than particle transition */
  double *lp0[sys->max_compound],*lp1[sys->max_compound];
  try{
    for(int c=0 ; c<sys->max_compound ; c++){
      lp0[c] = new double [MAX_LEVELS];
      lp1[c] = new double [MAX_LEVELS];
    }
  }
  catch(bad_alloc &e){
    message << "memory allocation error";
    cohTerminateCode("SpectraExclusive");
  }

  for(int c=0 ; c<sys->max_compound ; c++){
    for(int i=0; i< MAX_LEVELS ; i++) lp0[c][i] = lp1[c][i] = 0.0;
  }

  double renorm = specRenormReaction();

  for(int c0=0 ; c0<sys->max_compound ; c0++){
    int k0 = 0;
    if(c0 != 0) statStoreContinuumTransmission(c0,ncl[c0].excitation[k0],pdt,tc);

    spc->memclear("cn");

    for(k0=0 ; k0<ncl[c0].ncont ; k0++){
      statStoreContinuumGammaTransmission(k0,tg,&ncl[c0]);

      if((c0 == 0) && (k0 == 0) && (sys->incident.pid != gammaray)){
        specCompoundElastic(sys,renorm,tin,tc,td,tg,dir,spc);

        /*** Extract primary gamma-ray from level population */
        spc->npg = specStorePrimaryGamma(&ncl[c0],spc->pg,&gml);

        /*** For binary reactions, save discrete level population
             by particle only, and later subtract from the total population */
        for(int j=1 ; j<MAX_CHANNEL ; j++){
          if (!ncl[c0].cdt[j].status) continue;
          int p = ncl[c0].cdt[j].next;
          for(int i=0 ; i<ncl[p].ndisc ; i++) lp0[p][i] = ncl[p].lpop[i];
        }

      }else{
        statStoreDiscreteTransmission(c0,ncl[c0].excitation[k0],pdt,td);
        specCompoundDecay(c0,k0,tc,td,tg,spc);
      }

      /*** Population increment output */
      eclGenerateDecayTable(c0,k0,sys->max_channel,spc->dp);

      /*** Clear spectrum array */
      spc->memclear("dp");
    }

    /*** Particle decay from levels if possible */
    specLevelDecay(c0,pdt,&ncl[c0],tc,td,spc);

    /*** Save level population before cascading */
    for(int i=0 ; i<ncl[c0].ndisc ; i++) lp1[c0][i] = ncl[c0].lpop[i] - lp0[c0][i];

    /*** Remapping discrete transitions on to continuum bins */
    for(k0=ncl[c0].ncont ; k0<ncl[c0].ntotal ; k0++){
      specGammaCascadeExclusive(k0,spc->dp[0],&ncl[c0]);
      eclGenerateDecayTable(c0,k0,sys->max_channel,spc->dp);
      spc->memclear("dp");
    }
    specCumulativeSpectra(spc->getCsize(),spc->getNsize(),spc->cn,&ncl[c0]);
    crx.prod[c0].xsec += ncl[c0].lpop[0];
  }

  double reac = crx.total - crx.elastic;
  if(crx.total == 0.0) reac = crx.reaction; // for non-neutron reaction
  for(int c0=0 ; c0<sys->max_compound ; c0++){
    for(int i=0 ; i<ncl[c0].ndisc ; i++){
      lp0[c0][i] /= reac;
      lp1[c0][i] /= reac;
    }
  }
  eclLevelPopulation(sys->max_compound,lp0,lp1);

  for(int j=0 ; j<sys->max_compound ; j++){
    delete [] lp0[j];
    delete [] lp1[j];
  }
}


/**********************************************************/
/*      Renormalize Total Reaction Cross Section          */
/*      -----                                             */
/*             Subtract DWBA, PE, and DSD cross sections  */
/*             from the total reaction.                   */
/*             The CC cross sections, which are given by  */
/*             negative numbers, are not included         */
/**********************************************************/
double  specRenormReaction()
{
  double renorm = 1.0;

  double xs = 0.0;
  if(ctl.exciton) xs += crx.preeq;
  if(ctl.dsd)     xs += crx.dsd;
  for(int i=1 ; i<MAX_DIRECT ; i++){
    if(crx.direct[i] > 0.0) xs += crx.direct[i];
  }
  if(xs == 0.0) return(renorm);

  renorm = (crx.reaction-xs)/crx.reaction;
  if(renorm < 0.0){
    message << "the sum of direct, preeq etc. " << xs << " exceed total reaction " << crx.reaction;
    cohTerminateCode("specRenormReaction");
  }

  /*** Renormalize initial population calculated in statStoreInitPopulation */
  for(int j=0 ; j<MAX_J ; j++){
    ncl[0].pop[0][j].even *= renorm;
    ncl[0].pop[0][j].odd  *= renorm;
  }

  return(renorm);
}


/**********************************************************/
/*      For Monitoring Only                               */
/**********************************************************/
#ifdef DEBUG_POPCHECK
void specPopCheck(Nucleus *n)
{
  double s1 = 0.0, s2 = 0.0;

  cout << "#  ";
  cout << setw(3) << setfill('0') << n->za.getZ() << '-';
  cout << setw(3) << setfill('0') << n->za.getA() << endl;
  cout << setfill(' ');

  for(int k=0 ; k<n->ncont ; k++){
    cout.setf(ios::fixed, ios::floatfield);
    cout << setw(7) << setprecision(3) <<  n->excitation[k];
    cout.setf(ios::scientific, ios::floatfield);

    double s3 = 0.0;
    for(int j=0 ; j<10 ; j++) {
      s1 += n->pop[k][j].even+n->pop[k][j].odd;
      s3 += n->pop[k][j].even+n->pop[k][j].odd;
      cout << setprecision(2) << setw(9) << n->pop[k][j].even;
    }
    cout << endl;

    cout << "       ";
    for(int j=0 ; j<10 ; j++) {
      cout << setprecision(2) << setw(9) << n->pop[k][j].odd;
    }
    cout << setprecision(4) << setw(11) << s3 << endl;
  }

  s2 = 0.0;
  for(int k=0;k<n->ndisc;k++){
    s2 += n->lpop[k];
    cout << setw(5) << k
         << setprecision(4) << setw(13) << n->lev[k].energy
         << setw(13) << n->lpop[k] << endl;
  }

  cout << "# SUMcont " << s1 << endl;
  cout << "# SUMdisc " << s2 << endl;
  cout << "# SUMall  " << s1+s2 << endl;
  cout << endl;
  cout << endl;
}
#endif
