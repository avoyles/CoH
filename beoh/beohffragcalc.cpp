/******************************************************************************/
/*  beohffragcalc.cpp                                                         */
/*        Main calculation part for fission fragment decay                    */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "beoh.h"
#include "beohoutput.h"
#include "nucleus.h"
#include "masstable.h"
#include "beta2.h"

#define CALC_MONITOR

static bool PRINT_INDIV_RESULT        = false;
static bool PRINT_INDIV_SPEC          = false;
static bool PRINT_INDIV_CHI           = false;

static void   FFPMultiChanceYield (const int, const int);
static int    FFPMultiChanceFindIndex (const int, const int, const int, char *);
static int    FFPMultiChanceFindIndexL (const int, const int);
static int    FFPMultiChanceFindIndexH (const int, const int);
static void   FFPSetExcitationEnergy (const int, BData *, BData *, FissionFragmentPair *, const double);


static void   FFPSplitCalc1 (const int, const int, System *, System *, BData *, BData *, FFragData *, Pdata *);
static void   FFPSplitCalc2 (const int, System *, System *, BData *, BData *, Pdata *, const double);
static void   FFPPost (const int, System *, BData *, double, FragmentObservable *);
static void   FFPCombined (void);
static void   FFPCombinedSpectra (void);

#ifdef CALC_MONITOR
static void   FFPCalcMonitor (const int, const int, const int, const int, BData *, BData *);
#endif

extern FissionFragmentPair *ffp;
extern FissionObservable   fob;

static double yieldL = 0.0;
static double yieldH = 0.0;

/**********************************************************/
/*      For Each Fission Fragment Pair                    */
/**********************************************************/
void FFPSplit(CalcMode mode, System *sys, FFragData *fdt, const int ompid)
{
  BData       *bdtL, *bdtH;
  System      sysL, sysH;
  Pdata       pdt[MAX_CHANNEL];

  setupInitSystem(MAX_CHANNEL,pdt);
  pdt[neutron].omp = ompid;

  /*** Energy etc. data for both light and heavy fragments including multi-chance fission component */
  bdtL = new BData [fdt->getFissionChance()];
  bdtH = new BData [fdt->getFissionChance()];

  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    bdtL[nc].init();
    bdtH[nc].init();
    bdtL[nc].setMode(statdecay);
    bdtH[nc].setMode(statdecay);
    bdtL[nc].spinfactor = bdtH[nc].spinfactor = fdt->mc[nc].spinfactor;
    bdtL[nc].fraction   = bdtH[nc].fraction   = fdt->mc[nc].fraction;
  }

  /*** count number of total fragment pairs to be included */
  int ntotal = 0;
  for(int k=0 ; k<ffp[0].getN() ; k++){
    FFPMultiChanceYield(fdt->getFissionChance(),k);
    if( (yieldL + yieldH)*0.5 > fdt->ycutoff ) ntotal ++;
  }

  /*** loop over fission fragment pairs */
  int ncount = 1;
  for(int k=0 ; k<ffp[0].getN() ; k++){

//    if(ffp[0].fragment[k].getZl() != 37 || ffp[0].fragment[k].getAl() != 98 ) continue; // adhoc
//    if(ffp[0].fragment[k].getZh() != 53 || ffp[0].fragment[k].getAh() != 141 ) continue; // adhoc
//    if(ffp[0].fragment[k].getZh() != 52 || ffp[0].fragment[k].getAh() != 140 ) continue; // adhoc
//    if(ffp[0].fragment[k].getZh() != 52) continue;  // adhoc
//    if(ffp[0].fragment[k].getAh() != 138) continue;  // adhoc

    FFPMultiChanceYield(fdt->getFissionChance(),k);
    if( (yieldL + yieldH)*0.5 <= fdt->ycutoff ) continue;


    /*** Z,A for each fragment */
    sysL.init();
    sysH.init();

    sysL.compound.setZA(ffp[0].fragment[k].getZl(),ffp[0].fragment[k].getAl());
    sysH.compound.setZA(ffp[0].fragment[k].getZh(),ffp[0].fragment[k].getAh());

    /*** find beta2 deformation parameter */
    sysL.beta2 = beta2Read(&sysL.compound);
    sysH.beta2 = beta2Read(&sysH.compound);

    /*** copy energy bin width */
    sysL.energy_bin_in = sys->energy_bin_in;
    sysH.energy_bin_in = sys->energy_bin_in;

    /*** main calculation for each pair */
    if(mode != fissionspec){
      FFPSplitCalc1(k,fdt->getFissionChance(),&sysL,&sysH,bdtL,bdtH,fdt,pdt);
      /*** postprocess for L and H combined quantities */
      FFPCombined();
    }
    else{
      FFPSplitCalc2(k,&sysL,&sysH,&bdtL[0],&bdtH[0],pdt,fdt->mc[0].rt);
      /*** total spectra from L and H */
 /*
      cout.setf(ios::scientific, ios::floatfield);
      cout << "#  ";
      cout << setw(4) << ncount;
      cout << setw(4) << ffp[0].fragment[k].getZl();
      cout << setw(4) << ffp[0].fragment[k].getAl();
      cout << setprecision(5) << setw(14) << yieldL;
      cout << setw(4) << ffp[0].fragment[k].getZh();
      cout << setw(4) << ffp[0].fragment[k].getAh();
      cout << setprecision(5) << setw(14) << yieldH << endl;

      double de = (sys->energy_bin_in == 0.0) ? ENERGY_BIN : sys->energy_bin_in;

      for(int i=0 ; i<fob.getNsize() ; i++){
        double e = (i > 0)  ? (i - 0.5)*de : 0;
        cout << setprecision(5);
        cout << setw(14) << e;
        cout << setw(14) << fob.L.speclab[i];
        cout << setw(14) << fob.H.speclab[i] << endl;
        if(e > 30.0) break;
      }
      cout << endl;
      cout << endl;
*/
      FFPCombinedSpectra();
    }

#ifdef CALC_MONITOR
    /*** calculation monitor */
    FFPCalcMonitor(ntotal,ncount,k,fdt->getFissionChance(),bdtL,bdtH);
#endif

    /*** print individual result if needed */
    if(ncl[0].de > 0.0){
      if(PRINT_INDIV_SPEC) FFPOutputIndividualSpectrum(&ffp[0].fragment[k],fob.getNsize(),ncl[0].de,fob.L.spectrum[neutron],fob.H.spectrum[neutron],fob.spectrum[neutron]);

      if(PRINT_INDIV_CHI) FFPOutputIndividualSpectrum(&ffp[0].fragment[k],fob.getNsize(),ncl[0].de,fob.L.speclab,fob.H.speclab,fob.chi);
    }

    ncount ++;
  }

  /*** add prefission neutrons */
  fob.eprefis = 0.0;
  fob.nprefis = 0.0;
  for(int nc=1 ; nc < fdt->getFissionChance() ; nc++){
    fob.nprefis += nc * fdt->mc[nc].fraction;
    fob.eprefis += nc * fdt->mc[nc].fraction * fdt->mc[nc].eprefis;
  }


  delete [] bdtL;
  delete [] bdtH;
}


/**********************************************************/
/*      Fragment Yield Including MultiChance Components   */
/**********************************************************/
void FFPMultiChanceYield(const int nc, const int k0)
{
  if(nc == 1){
    yieldL = yieldH = ffp[0].fragment[k0].yield;
    return;
  }
  else{
    yieldL = yieldH = ffp[0].fragment[k0].yield * ffp[0].fraction;
  }

  char z;
  for(int i=1 ; i<nc ; i++){
    /*** for the light fragemnt */
    int kl = FFPMultiChanceFindIndex(i,ffp[0].fragment[k0].getZl(),ffp[0].fragment[k0].getAl(),&z);
    if(kl >= 0) yieldL += ffp[i].fragment[kl].yield * ffp[i].fraction;
    /*** for the heavy fragemnt */
    int kh = FFPMultiChanceFindIndex(i,ffp[0].fragment[k0].getZh(),ffp[0].fragment[k0].getAh(),&z);
    if(kh >= 0) yieldH += ffp[i].fragment[kh].yield * ffp[i].fraction;
  }
}


/**********************************************************/
/*      Index for Same Z,A in Multi-Chance List           */
/**********************************************************/
int FFPMultiChanceFindIndex(const int c, const int z0, const int a0, char *d)
{
  ZAnumber za0(z0,a0);

  int kx = -1;
  *d = 'X';

  for(int k1=0 ; k1<ffp[c].getN() ; k1++){
    ZAnumber zaL(ffp[c].fragment[k1].getZl(),ffp[c].fragment[k1].getAl());
    ZAnumber zaH(ffp[c].fragment[k1].getZh(),ffp[c].fragment[k1].getAh());

    if(za0 == zaL){ kx = k1; *d = 'L'; break;}
    if(za0 == zaH){ kx = k1; *d = 'H'; break;}
  }
  return kx;
}


int FFPMultiChanceFindIndexL(const int c, const int k0)
{
  int k = -1;
  for(int k1=0 ; k1<ffp[c].getN() ; k1++){
    if( (ffp[0].fragment[k0].getZl() == ffp[c].fragment[k1].getZl()) &&
        (ffp[0].fragment[k0].getAl() == ffp[c].fragment[k1].getAl()) ){ k = k1; break; }
  }
  return k;
}

int FFPMultiChanceFindIndexH(const int c, const int k0)
{
  int k = -1;
  for(int k1=0 ; k1<ffp[c].getN() ; k1++){
    if( (ffp[0].fragment[k0].getZh() == ffp[c].fragment[k1].getZh()) &&
        (ffp[0].fragment[k0].getAh() == ffp[c].fragment[k1].getAh()) ){ k = k1; break; }
  }
  return k;
}


/**********************************************************/
/*      Fragment Pair for Yield Calculation               */
/*      Energy integration performed by energy            */
/*      distribution of initial population                */
/**********************************************************/
void FFPSplitCalc1(const int k, const int nc, System *sysL, System *sysH, BData *bdtL, BData *bdtH, FFragData *fdt, Pdata *pdt)
{
  BData bL, bH;
  /*** for the first chance fission */
  FFPSetExcitationEnergy(k,&bL,&bH,&ffp[0],fdt->mc[0].rt);
  bdtL[0].exmean  = sysL->ex_total = bL.exmean;
  bdtH[0].exmean  = sysH->ex_total = bH.exmean;
  bdtL[0].exwidth = bL.exwidth;
  bdtH[0].exwidth = bH.exwidth;
  bdtL[0].ekmean  = bL.ekmean;
  bdtH[0].ekmean  = bH.ekmean;

  /*** multi-chance fission case */
  char z;
  for(int nc=1 ; nc < fdt->getFissionChance() ; nc++){

    /*** find the same ZA in the multi-chance chain */
    int kl = FFPMultiChanceFindIndex(nc,ffp[0].fragment[k].getZl(),ffp[0].fragment[k].getAl(),&z);
    FFPSetExcitationEnergy(kl,&bL,&bH,&ffp[nc],fdt->mc[nc].rt);
    if(z == 'L'){
      bdtL[nc].exmean  = bL.exmean;  if(bL.exmean > sysL->ex_total) sysL->ex_total = bL.exmean;
      bdtL[nc].exwidth = bL.exwidth;
      bdtL[nc].ekmean  = bL.ekmean;
    }
    else{
      bdtL[nc].exmean  = bH.exmean;  if(bH.exmean > sysL->ex_total) sysL->ex_total = bH.exmean;
      bdtL[nc].exwidth = bH.exwidth;
      bdtL[nc].ekmean  = bH.ekmean;
    }

    int kh = FFPMultiChanceFindIndex(nc,ffp[0].fragment[k].getZh(),ffp[0].fragment[k].getAh(),&z);
    FFPSetExcitationEnergy(kh,&bL,&bH,&ffp[nc],fdt->mc[nc].rt);
    if(z == 'L'){
      bdtH[nc].exmean  = bL.exmean;  if(bL.exmean > sysH->ex_total) sysH->ex_total = bL.exmean;
      bdtH[nc].exwidth = bL.exwidth;
      bdtH[nc].ekmean  = bL.ekmean;
    }
    else{
      bdtH[nc].exmean  = bH.exmean;  if(bH.exmean > sysH->ex_total) sysH->ex_total = bH.exmean;
      bdtH[nc].exwidth = bH.exwidth;
      bdtH[nc].ekmean  = bH.ekmean;
    }
  }

  /*** exec Hauser-Feshbach for each fragment, post processing, collecting calculated results */
  fob.L.clear();
  beohFFDecay(nc,sysL,pdt,bdtL,PRINT_INDIV_RESULT);
  FFPPost(nc,sysL,bdtL,yieldL,&fob.L);

  fob.H.clear();
  beohFFDecay(nc,sysH,pdt,bdtH,PRINT_INDIV_RESULT);
  FFPPost(nc,sysH,bdtH,yieldH,&fob.H);
}


/**********************************************************/
/*      Store Average Excitation Energies                 */
/**********************************************************/
void FFPSetExcitationEnergy(const int k, BData *bL, BData *bH, FissionFragmentPair *f, const double rt)
{
  if(k < 0){
    bL->exmean = bL->exwidth = bL->ekmean = 0.0;
    bH->exmean = bH->exwidth = bH->ekmean = 0.0;
    return;
  }

  /*** for the first chance fission */
  ZAnumber cL, cH;
  cL.setZA(f->fragment[k].getZl(),f->fragment[k].getAl());
  cH.setZA(f->fragment[k].getZh(),f->fragment[k].getAh());

  /*** energy sorting by the RT model */
  double eratio = beohExcitationEnergyDivide(&cL,&cH,rt,f->fragment[k].ex);

  if(eratio > 0.0){
    bL->exmean = f->fragment[k].ex * eratio;
    bH->exmean = f->fragment[k].ex * (1.0 - eratio);
  }
  else{
    bL->exmean = 0.0;
    bH->exmean = 0.0;
  }

  /*** when dE is given for TXE, divide dE into two fragments */
  if(f->fragment[k].sigma > 0.0){
    double c = ffp[0].fragment[k].sigma / sqrt(bL->exmean * bL->exmean + bH->exmean * bH->exmean);
    bL->exwidth = c * bL->exmean;
    bH->exwidth = c * bH->exmean;
  }
  else{
    bL->exwidth = 0.0;
    bH->exwidth = 0.0;
  }

  /*** average kinetic energies */
  double x = f->fragment[k].getAl() + f->fragment[k].getAh();
  bL->ekmean = f->fragment[k].ek * f->fragment[k].getAh() / x;
  bH->ekmean = f->fragment[k].ek * f->fragment[k].getAl() / x;
}


/**********************************************************/
/*      Fragment Pair for Fission Spectrum Calculation    */
/*      Perform exact energy intergration                 */
/**********************************************************/
void FFPSplitCalc2(const int k, System *sysL, System *sysH, BData *bdtL, BData *bdtH, Pdata *pdt, double rt)
{
  double *speclab = new double [fob.getNsize()];
  double *tkedist = new double [fob.getKsize()];

  for(int j=0 ; j<fob.getNsize() ; j++){
    fob.L.speclab[j] = fob.H.speclab[j] = speclab[j] = 0.0;
    for(int p=0 ; p<2 ; p++) fob.L.spectrum[p][j] = fob.H.spectrum[p][j] = 0.0;
  }

  for(int i=0 ; i<fob.getKsize() ; i++){
    fob.L.nTKE[i] = fob.L.yTKE[i] = fob.H.nTKE[i] = fob.H.yTKE[i] = tkedist[i] = 0.0;
  }

  /*** determine Jmax for each fragment at the mean energy */
  FFPSetExcitationEnergy(k,bdtL,bdtH,&ffp[0],rt);
  sysL->ex_total = bdtL->exmean;
  sysH->ex_total = bdtH->exmean;

  int jmaxL = beohFFDecayPrep(sysL,pdt);
  int jmaxH = beohFFDecayPrep(sysH,pdt);

  /*** kinetic energy increment */
  int kmax = ffp[0].fragment[k].ek + EXRANGE_FACTOR * ffp[0].fragment[k].sigma;
  int kmin = ffp[0].fragment[k].ek - EXRANGE_FACTOR * ffp[0].fragment[k].sigma; if(kmin < 0) kmin = 0;
  int nk   = (double)(kmax - kmin) + 1;

  /*** Gaussian weight */
  double sg = 0.0;
  for(int i=0 ; i<nk ; i++){
    double ek = kmin + i;
    if(ffp[0].fragment[k].ex + ffp[0].fragment[k].ek - ek < 0.0) {nk = i; break;}
    tkedist[i] = exp(-(ek - ffp[0].fragment[k].ek)*(ek - ffp[0].fragment[k].ek)
               /(2.0*ffp[0].fragment[k].sigma*ffp[0].fragment[k].sigma));
    sg += tkedist[i];
  }
  if(sg > 0.0){ for(int i=0 ; i<nk ; i++) tkedist[i] /= sg; }

  /*** for each TKE bin */
  double nubarL = 0.0, nubarH = 0.0;
  for(int i=0 ; i<nk ; i++){
    double tke = (double)(kmin + i);
    double txe = ffp[0].fragment[k].ex + ffp[0].fragment[k].ek - tke;

//  if(tke < 100.0 || tke > 120.0) continue; // adhoc

    /*** divide TKE and TXE into two fragments */
    double eratio = beohExcitationEnergyDivide(&sysL->compound,&sysH->compound,rt,txe);
    bdtL->exmean = txe * eratio;
    bdtH->exmean = txe * (1.0 - eratio);

    bdtL->ekmean = tke * sysH->compound.getA() / (sysL->compound.getA() + sysH->compound.getA());
    bdtH->ekmean = tke * sysL->compound.getA() / (sysL->compound.getA() + sysH->compound.getA());

    /*** statistical decay calculation and postprocessing for the light fragment */
    beohFFDecaySpec(sysL,pdt,bdtL,jmaxL,fob.getNsize(),speclab);
    beohAverageEnergy(MAX_ENERGY_BIN,fob.L.multiplicity,fob.L.eaverage,fob.L.etotal,crx.spectra);
    for(int j=0 ; j<fob.getNsize() ; j++){
      for(int p=0 ; p<2 ; p++) fob.L.spectrum[p][j] += crx.spectra[p][j] * tkedist[i];
      fob.L.speclab[j] += speclab[j] * tkedist[i];
    }

    /*** statistical decay calculation and postprocessing for the heavy fragment */
    beohFFDecaySpec(sysH,pdt,bdtH,jmaxH,fob.getNsize(),speclab);
    beohAverageEnergy(MAX_ENERGY_BIN,fob.H.multiplicity,fob.H.eaverage,fob.H.etotal,crx.spectra);
    for(int j=0 ; j<fob.getNsize() ; j++){
      for(int p=0 ; p<2 ; p++) fob.H.spectrum[p][j] += crx.spectra[p][j] * tkedist[i];
      fob.H.speclab[j] += speclab[j] * tkedist[i];
    }

    fob.L.nTKE[kmin+i] = fob.L.multiplicity[neutron];
    fob.H.nTKE[kmin+i] = fob.H.multiplicity[neutron];

    nubarL += fob.L.multiplicity[neutron] * tkedist[i];
    nubarH += fob.H.multiplicity[neutron] * tkedist[i];
  }

  fob.L.multiplicity[neutron] = nubarL;
  fob.H.multiplicity[neutron] = nubarH;

  for(int i=0 ; i<nk ; i++){
    fob.yTKE[kmin+i] += ffp[0].fragment[k].yield;
    fob.nTKE[kmin+i] += (fob.L.nTKE[kmin+i] + fob.H.nTKE[kmin+i]) * ffp[0].fragment[k].yield;
  }

  /*** just in case if TKE look above was skipped */
  cohResetNucleus();
  cohResetNuclearStructure();

  delete [] speclab;
  delete [] tkedist;
}



/**********************************************************/
/*      Accumulate Decay Calculation Results              */
/**********************************************************/
void FFPPost(const int nc, System *sys, BData *bdt, const double yield, FragmentObservable *fo)
{
  /*** photon and neutron spectra */
  for(int p=0 ; p<2 ; p++){
    for(int i=0 ; i<fob.getNsize() ; i++) fo->spectrum[p][i] = crx.spectra[p][i];
  }

  /*** simple conversion of neutron cms-spectrum to lab */
  for(int i=0 ; i<fob.getNsize() ; i++) fo->speclab[i] = crx.spectra[neutron][i];
  beohIntegCmsSpectrum(fob.getNsize(),nc,sys->compound.getA(),bdt,ncl[0].de,fo->speclab);

  /*** average energy and multiplicity */
  beohAverageEnergy(fob.getNsize(),fo->multiplicity,fo->eaverage,fo->etotal,fo->spectrum);

  /*** neutron multiplicity distributions, P(nu) */
  for(int i=0 ; i<fob.getLsize() ; i++) fo->Pn[i] = 0.0;
  int a0 = ncl[0].za.getA();
  for(int j=0 ; j<sys->max_compound ; j++){
    int da = a0 - ncl[j].za.getA();
    if(da < fob.getLsize()) fo->Pn[da] = yield * crx.prod[j].xsec;
  }

  /*** multiplicity and average energy distributions as function of A */
  for(int p=0 ; p<2 ; p++){
    fob.multiplicitydist[p][sys->compound.getA()] += yield * fo->multiplicity[p];
    fob.eaveragedist[p][sys->compound.getA()]     += yield * fo->eaverage[p];
  }

  /*** pre neutron emission chain yield */
  fob.preChainYield[sys->compound.getA()] += yield;

  /*** post neutron emission chain yield */
  for(int j=0 ; j<sys->max_compound ; j++){
    fob.postChainYield[ncl[j].za.getA()] += yield * crx.prod[j].xsec;
  }

  /*** add long-lived isotope productions */
  outPrepCumulativeResidual(sys->max_compound,yield);
}


/**********************************************************/
/*      Postprocess for L and H Combined Quantities       */
/**********************************************************/
void FFPCombined(void)
{
  /*** average photon and neutron multiplicities */
  for(int p=0 ; p<2 ; p++){
    double mtot = fob.L.multiplicity[p] + fob.H.multiplicity[p];
    double etot = yieldL * fob.L.multiplicity[p] * fob.L.eaverage[p] + yieldH * fob.H.multiplicity[p] * fob.H.eaverage[p];
    double eave = (mtot > 0.0) ? etot / mtot : 0.0;

    fob.multiplicity[p] += yieldL * fob.L.multiplicity[p] + yieldH * fob.H.multiplicity[p];
    fob.eaverage[p]     += eave;
    fob.etotal[p]       += etot;
  }

  /*** calculate neutron emission probability distribution
       P(0) = Pl(0) * Ph(0), P(1) = Pl(0) * Ph(1) + Pl(1) * Ph(0),
       P(2) = Pl(0) * Ph(2) + Pl(1) * Ph(1) + Pl(0) * Ph(2), ... */
  for(int i=0 ; i<fob.getLsize() ; i++){
    for(int j=0 ; j<=i ; j++){
      fob.Pn[i] += fob.L.Pn[j] * fob.H.Pn[i-j];
    }
  }

  /*** average spectra */
  FFPCombinedSpectra();
}


void FFPCombinedSpectra(void)
{
  /*** average spectra */
  for(int i=0 ; i<fob.getNsize() ; i++){
    for(int p=0 ; p<2 ; p++){
      fob.L.spectrum[p][i] *= yieldL;
      fob.H.spectrum[p][i] *= yieldH;
      fob.spectrum[p][i]   += fob.L.spectrum[p][i] + fob.H.spectrum[p][i];
    }
    fob.L.speclab[i] *= yieldL;
    fob.H.speclab[i] *= yieldH;
    fob.chi[i]       += fob.L.speclab[i] + fob.H.speclab[i];
  }
}


/**********************************************************/
/*      Calculation Monitor                               */
/**********************************************************/
#ifdef CALC_MONITOR
static bool first_call= true;
void FFPCalcMonitor(const int m, const int n, const int k, const int nc, BData *bl, BData *bh)
{
  static bool detailed = true;

  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  if(first_call){
    cout << "# Ntotal = " << m << endl;
    cout << "#    n    k  Zl  Al  Zh  Ah";
    cout << "  YieldL      YieldH      TKE[MeV]    TXE[MeV]    El[MeV]     Wl[MeV]     Eh[MeV]     Wh[MeV]" << endl;
    first_call = false;
  }

  double ek = 0.0, ex = 0.0;
  double blm = 0.0, blw = 0.0, bhm = 0.0, bhw = 0.0;
  if(nc == 1){
    ex = ffp[0].fragment[k].ex;
    ek = ffp[0].fragment[k].ek;
    blm = bl[0].exmean;
    bhm = bh[0].exmean;
    blw = bl[0].exwidth;
    bhw = bh[0].exwidth;
  }
  else{
    for(int i=0 ; i<nc ; i++){
      ex  += ffp[i].fraction * ffp[i].fragment[k].ex;
      ek  += ffp[i].fraction * ffp[i].fragment[k].ek;
      blm += ffp[i].fraction * bl[i].exmean;
      bhm += ffp[i].fraction * bh[i].exmean;
      blw += ffp[i].fraction * bl[i].exwidth * bl[i].exwidth;
      bhw += ffp[i].fraction * bh[i].exwidth * bh[i].exwidth;
    }
    blw = sqrt(blw);
    bhw = sqrt(bhw);
  }

  cout << "# " << setw(4) << n << setw(5) << k;
  cout << setw(4) << ffp[0].fragment[k].getZl() << setw(4) << ffp[0].fragment[k].getAl();
  cout << setw(4) << ffp[0].fragment[k].getZh() << setw(4) << ffp[0].fragment[k].getAh();
  cout << setw(12) << yieldL;
  cout << setw(12) << yieldH;
  cout << setw(12) << ek;
  cout << setw(12) << ex;
  cout << setw(12) << blm;
  cout << setw(12) << blw;
  cout << setw(12) << bhm;
  cout << setw(12) << bhw << endl;

  if(nc > 1 && detailed){
    for(int i=0 ; i<nc ; i++){
      cout << "#              ";
      cout << setw(12) << ffp[i].fraction;

      int kl = FFPMultiChanceFindIndexL(i,k);
      int kh = FFPMultiChanceFindIndexH(i,k);
      cout << setw(12) << ffp[i].fragment[kl].yield;
      cout << setw(12) << ffp[i].fragment[kh].yield;
      cout << setw(12) << ffp[i].fragment[k].ek;
      cout << setw(12) << ffp[i].fragment[k].ex;
      cout << setw(12) << bl[i].exmean;
      cout << setw(12) << bl[i].exwidth;
      cout << setw(12) << bh[i].exmean;
      cout << setw(12) << bh[i].exwidth << endl;
    }
  }
}
#endif

