/******************************************************************************/
/*  beohffrag.cpp                                                             */
/*        Calculation of statistical decay of fission fragments               */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "beoh.h"
#include "beohoutput.h"
#include "ffpmultichance.h"
#include "nucleus.h"
#include "masstable.h"
#include "terminate.h"

static bool PRINT_FISSION_PAIR        = false;
static bool PRINT_FISSION_PAIR_SUPPL  = false;

static void beohFissionFragmentInit (System *, BData *, FFragData *);
static void beohEnergyDependentParameters (FFragData *);
static void beohFFPAllocateMemory (const int);
static void beohFFPDeleteAllocated (const int);


/**********************************************************/
/*      Global Scope Data Arrays                          */
/**********************************************************/
/* data for fission yield Y(Z,A,TXE) */
FissionFragmentPair *ffp;

/* data for post-processing */
FissionObservable fob;        // post-processed quantities for data comparison


/**********************************************************/
/*      Fission Fragment Decay Calculation                */
/**********************************************************/
void beohFissionFragmentMain(System *sys, BData *bdt, FFragData *fdt, unsigned int ompid)
{
  if(fdt->getFissionChance() < 1) return;

  if( (bdt->getMode() == fissionspec) && (fdt->getFissionChance() > 1) ){
    message << "fission specturm calculation for multi-chance case not supported";
    cohTerminateCode("beohFissionFragmentMain");
  }

  if(fdt->ycutoff == 0.0) fdt->ycutoff = YIELD_CUTOFF;
  if(ompid == 0) ompid = 6; // Koning-Delaroche potential as default


  /*** allocate Fragment Pair Object for each chance-fission */
  ffp = new FissionFragmentPair [fdt->getMaxFissionChance()];

  /*** determine fissioning system for the first chance fission and set binding energies */
  beohFissionFragmentInit(sys,bdt,fdt);
  beohEnergyDependentParameters(fdt);

  /*** prepare fission fragment pair object and allocate memory */
  beohFFPAllocateMemory(fdt->getMaxFissionChance());
  beohFFDecayAllocateMemory();

  /*** create fission fragment pairs and store them in the FissionFragmentPair object */
  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){

    /*** 'ffydata'== input, Wahl or internal */
    if(fdt->mc[nc].ffydata == "input" || fdt->mc[nc].ffydata == "Wahl" || fdt->mc[nc].ffydata == "internal"){

      /*** retrieve Gaussian parameters*/
      FFPMassChainYieldParameters(fdt->spontaneous,fdt->mc[nc].ffydata,fdt->mc[nc].GaussSigma,fdt->mc[nc].GaussDelta,fdt->mc[nc].GaussFract,&ffp[nc]);

      /*** generate fission fragment pairs and re-calculate TKE */
      FFPGeneratePairs(fdt->spontaneous,fdt->mc[nc].tke,&ffp[nc],fdt->mc[nc].ZpFactor);
    }

    /*** Raw data file method */
    else{
      /*** read data in external file */
      FFPReadRawYield(fdt->mc[nc].ffydata,&ffp[nc]);
    }

    /*** recalculate TKE from fission product pairs, only for printing */
    fdt->mc[nc].tke = FFPAverageTKE(&ffp[nc]);
  }

  /*** ffp[0] will have a union list of the all FF appearing at all multi-chances */
  if(fdt->getFissionChance() > 1) FFPExpandList(fdt->getFissionChance(),ffp);

  /*** calculate decay of fission fragments */
  FFPSplit(bdt->getMode(),sys,fdt,ompid);

  /*** print parameters */
  if(PRINT_FISSION_PAIR){
    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++) FFPOutputPairs(fdt->ycutoff,&ffp[nc]);
  }
  outSystem(bdt->getMode(),sys);
  outFissionFragment(sys->compound.getA(),fdt);

  /*** main output */
  double de = (sys->energy_bin_in == 0.0) ? ENERGY_BIN : sys->energy_bin_in;

  if(bdt->getMode() == fissionspec) FFPOutputSpec(de,fdt->maxtemperature,&fob);
  else{
    FFPOutput(de,fdt->massresolution,fdt->maxtemperature,&fob);

    /*** print production cross sections of ground and isomeric states */
    if(bdt->getMode() == cumulativeyield){
      beohFragmentBetaDecay(fdt->maxhalflife,fdt->branchdatafile);
    }
    else{
      outPrepPrintResidual();
      outPrepPrintIsomericRatio();
    }
  }

  /*** supplemental output */
  if(PRINT_FISSION_PAIR_SUPPL){
    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
      FFPOutputChainYield(ffp[0].nmass,ffp[nc].ac,ffp[nc].mass_first,fdt->mc[nc].GaussSigma,fdt->mc[nc].GaussDelta,fdt->mc[nc].GaussFract);
      FFPOutputTKEDistribution(3,&ffp[nc]);
      FFPOutputTXEDistribution(&ffp[nc]);
    }
  }

  /*** free allocated memory */
  beohFFPDeleteAllocated(fdt->getMaxFissionChance());
  beohFFDecayDeleteAllocated();

  delete [] ffp;
}


/**********************************************************/
/*      Determine Fissioning System                       */
/**********************************************************/
void beohFissionFragmentInit(System *sys, BData *bdt, FFragData *fdt)
{
  /*** neutron induced fission case */
  if(bdt->energy >= 0.0){
    ffp[0].en = bdt->energy;

    /*** first-chance fissioning nucleus = A+1 */
    sys->compound.setZA(sys->target.getZ(),sys->target.getA()+1);

    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
      /*** ZA of fissioning compound nucleus */
      ffp[nc].setZA(sys->compound.getZ(),sys->compound.getA() - nc);

      /*** neutron binding energy */
      ffp[nc].bn = mass_excess(ffp[nc].zf,ffp[nc].af-1) + ENEUTRON
                 - mass_excess(ffp[nc].zf,ffp[nc].af  );
    }

    /*** total excitation energy Bn + En(cms) */
    sys->ex_total = ffp[0].en * (double)sys->target.getA()/(sys->target.getA()+1.0) + ffp[0].bn;

    /*** multi-chance fission data from file */
    if(fdt->mcffile != ""){
      double dat[16];
      int ncdat = FFPReadMultichanceFissionData(fdt->mcffile,ffp[0].en,dat);
      for(int nc=0 ; nc < ncdat ; nc++){
        if(fdt->mc[nc].fraction == 1.0) fdt->mc[nc].fraction = dat[nc    ];
        if(fdt->mc[nc].tke      == 0.0) fdt->mc[nc].tke      = dat[nc + 4];
        if(fdt->mc[nc].exfis    == 0.0) fdt->mc[nc].exfis    = dat[nc + 8];
        if(fdt->mc[nc].eprefis  == 0.0) fdt->mc[nc].eprefis  = dat[nc +12];
      }
      /*** when data file contains more fission-chance, copy the first chance data */
      if(fdt->getFissionChance() < ncdat){
        for(int nc=fdt->getFissionChance() ; nc < ncdat ; nc++){
          ffp[nc].setZA(sys->compound.getZ(),sys->compound.getA() - nc);
          ffp[nc].bn = mass_excess(ffp[nc].zf,ffp[nc].af-1) + ENEUTRON
                     - mass_excess(ffp[nc].zf,ffp[nc].af  );
          fdt->mc[nc].spinfactor = fdt->mc[0].spinfactor;
          fdt->mc[nc].rt         = fdt->mc[0].rt;
          fdt->mc[nc].ffydata    = fdt->mc[0].ffydata;
        }
      }
      fdt->resetFissionChance(ncdat);
    }

    /*** enforce single-chance calculation when Pf=1.0 */
    if(fdt->mc[0].fraction == 1.0) fdt->resetFissionChance(1);

    /*** fission probabilities */
    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
      ffp[nc].fraction = fdt->mc[nc].fraction;
      if(ffp[nc].fraction == 0.0){ fdt->resetFissionChance(nc); break; }
    }

    /*** check closed channel and average exication energy causing fission */
    double bnsum = 0.0;
    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){

      /*** when not given, we assume fission happens at the highest excitation energy */
      if(fdt->mc[nc].exfis == 0.0){ fdt->mc[nc].exfis = sys->ex_total - bnsum; }
      bnsum += ffp[nc].bn;

      ffp[nc].ex = fdt->mc[nc].exfis;
      /*** check closed channel */
      if(ffp[nc].ex < 0.0){ fdt->resetFissionChance(nc); break; }
    }

    fdt->spontaneous = false;
  }
  /*** spontaneous fission case */
  else{
    sys->compound = sys->target;
    ffp[0].en = ffp[0].ex = 0.0;

    fdt->spontaneous = true;
  }

  /*** check TKE */
  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    if(fdt->mc[nc].tke == 0.0){
      fdt->mc[nc].tke = FFPSystematics_TKE_En(fdt->spontaneous,ffp[nc].zf,ffp[nc].af,ffp[nc].ex,ffp[nc].bn);
    }
  }

  /*** determine mass-range to be calculated */
  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    ffp[nc].setZARange(MASS_FIRST,CHARGE_RANGE);
  }
}


/**********************************************************/
/*      Calculate Energy-Dependent Parameters             */
/**********************************************************/
void beohEnergyDependentParameters(FFragData *fdt)
{
  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    double en = ffp[nc].ex - ffp[nc].bn;
    if(en < 0.0) en = 0.0;

    if(fdt->mc[nc].tke1 != 0.0){
      fdt->mc[nc].tke = fdt->mc[nc].tke + en * fdt->mc[nc].tke1;
      if(fdt->mc[nc].tke < 0.0) fdt->mc[nc].tke = 0.0;
    }

    /*** set lower limit of fJ to be  0.1 */
    if(fdt->mc[nc].spinfactor1 != 0.0){
      fdt->mc[nc].spinfactor  = fdt->mc[nc].spinfactor + en * fdt->mc[nc].spinfactor1;
      if(fdt->mc[nc].spinfactor < 0.1) fdt->mc[nc].spinfactor = 0.1;
    }

    /*** when RT < 1, then RT is always 1.0 */
    if(fdt->mc[nc].rt1 != 0.0){
      fdt->mc[nc].rt  = fdt->mc[nc].rt + en * fdt->mc[nc].rt1;
      if(fdt->mc[nc].rt < 1.0) fdt->mc[nc].rt = 1.0;
    }

    /*** set Zp parameters lower limit 1e-10, since zero does not change calculations */
    for(int i=0 ; i<2 ; i++){
      if(fdt->mc[nc].ZpFactor[i+2] != 0.0){
        fdt->mc[nc].ZpFactor[i] = fdt->mc[nc].ZpFactor[i] + en * fdt->mc[nc].ZpFactor[i+2];
        if(fdt->mc[nc].ZpFactor[i] <= 0.0) fdt->mc[nc].ZpFactor[i] = 1e-10;
      }
    }

    for(int i=0 ; i<4 ; i++){
      if(fdt->mc[nc].GaussFract[i+4] != 0.0){
        fdt->mc[nc].GaussFract[i] = fdt->mc[nc].GaussFract[i] + en * fdt->mc[nc].GaussFract[i+4];
        if(fdt->mc[nc].GaussFract[i] < 0.0){
          message << i+1 << "-th primary gaussian fraction became negative";
          cohTerminateCode("beohEnergyDependentParameters");
        }
      }
      if(fdt->mc[nc].GaussSigma[i+4] != 0.0){
        fdt->mc[nc].GaussSigma[i] = fdt->mc[nc].GaussSigma[i] + en * fdt->mc[nc].GaussSigma[i+4];
        if(fdt->mc[nc].GaussSigma[i] < 0.0){
          message << i+1 << "-th primary gaussian width became negative";
          cohTerminateCode("beohEnergyDependentParameters");
        }
      }
      if(i == 3) continue; // symmetric part does not have mass shift
      if(fdt->mc[nc].GaussDelta[i+4] != 0.0){
        fdt->mc[nc].GaussDelta[i] = fdt->mc[nc].GaussDelta[i] + en * fdt->mc[nc].GaussDelta[i+4];
        if(fdt->mc[nc].GaussDelta[i] < 0.0){
          message << i+1 << "-th primary gaussian mass shift became negative";
          cohTerminateCode("beohEnergyDependentParameters");
        }
      }
    }
  }
}


/**********************************************************/
/*      Allocate Arrays                                   */
/**********************************************************/
void beohFFPAllocateMemory(const int nc)
{
  try{
    for(int i=0 ; i<nc ; i++) ffp[i].memalloc(MAX_FISSION_PAIR);
    fob.memalloc(MAX_MULTIPLICITY_DIST,ffp[0].mass_first+ffp[0].nmass+1,MAX_ENERGY_BIN);

    /*** allocate long-lived states */
    outPrepAllocateResidual(MAX_FISSION_PAIR * 2);
  }
  catch(bad_alloc &e){
    message << "memory allocation error";
    cohTerminateCode("beohFFPAllocateMemory");
  }
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void beohFFPDeleteAllocated(const int nc)
{
  for(int i=0 ; i<nc ; i++) ffp[i].memfree();
  fob.memfree();
  outPrepFreeResidual();
}

