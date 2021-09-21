/******************************************************************************/
/*  ffpoutput.cpp                                                             */
/*        Print Z, A, TKE, TXE, Y(A,Z)                                        */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "beoh.h"
#include "outformat.h"

static void FFPOutputMultiplicity (double, double *, double *, double *, double, double);
static void FFPOutputPnu (const int, double *);
static void FFPOutputYield (const int, double *, double *, double **, double **);
static void FFPOutputnuTKE (const int, double *);
static void FFPOutputSpectrum (const int, const int, const double, const double, double **, double *);
static void FFPFindARange (FissionFragmentPair *);
static void FFPGaussianBroadening (const int, const double, double *);

static int ahmin = 1000, ahmax = 0;

/**********************************************************/
/*      Output Calculated Results                         */
/**********************************************************/
void FFPOutput(const double de, const double dmass, const double maxtemp, FissionObservable *fob)
{
  /*** normalization by pre-neutron emission yield */
  /*** independent yields and quantities as a function of mass */
  double ytot = 0.0;
  for(int i=0 ; i<fob->getMsize() ; i++){
    if(fob->preChainYield[i] > 0.0){
      for(int p=0 ; p<2 ; p++){
        fob->multiplicitydist[p][i] /= fob->preChainYield[i];
        fob->eaveragedist[p][i]     /= fob->preChainYield[i];
      }
    }
    else{
      for(int p=0 ; p<2 ; p++){
        fob->multiplicitydist[p][i] = 0.0;
        fob->eaveragedist[p][i]     = 0.0;
      }
    }
    ytot += fob->preChainYield[i];
  }

  /*** when mass resolution is given, postY is broadened */
  if(dmass > 0.0){
    FFPGaussianBroadening(fob->getMsize(),dmass,fob->postChainYield);
    for(int p=0 ; p<2 ; p++){
      FFPGaussianBroadening(fob->getMsize(),dmass,fob->multiplicitydist[p]);
      FFPGaussianBroadening(fob->getMsize(),dmass,fob->eaveragedist[p]);
    }
  }

  /*** normalize Pnu */
  double x = 0.0;
  for(int i=0 ; i<fob->getLsize() ; i++) x += fob->Pn[i];
  if(x > 0.0){
    for(int i=0 ; i<fob->getLsize() ; i++) fob->Pn[i] /= x;
  }

  /*** normalize LAB neutron spectrum */
  x = 0.0;
  for(int k=0 ; k<fob->getNsize() ; k++) x += fob->chi[k] * de;
  if(x > 0.0){
    for(int k=0 ; k<fob->getNsize() ; k++) fob->chi[k] /= x;
  }


  /*** average multiplicity */
  FFPOutputMultiplicity(ytot,fob->multiplicity,fob->eaverage,fob->etotal,fob->eprefis,fob->nprefis);

  /*** print multiplicity distributions */
  FFPOutputPnu(fob->getLsize(),fob->Pn);

  /*** print yield, mulitplicity, and average energy distributions */
  FFPOutputYield(fob->getMsize(),fob->preChainYield,fob->postChainYield,fob->multiplicitydist,fob->eaveragedist);

  /*** print photon and neutron spectra */
  FFPOutputSpectrum(beohZeroCut(fob->spectrum[gammaray]),beohZeroCut(fob->chi),de,maxtemp,fob->spectrum,fob->chi);
}


/**********************************************************/
/*      Output TKE-dist Related Results and Chi           */
/**********************************************************/
void FFPOutputSpec(const double de, const double maxtemp, FissionObservable *fob)
{
  /*** normalize nu(TKE) */
  for(int i=0 ; i<fob->getKsize() ; i++){
    if(fob->yTKE[i] > 0.0) fob->nTKE[i] /= fob->yTKE[i];
  }

  /*** normalize LAB neutron spectrum */
  double x = 0.0;
  for(int k=0 ; k<fob->getNsize() ; k++) x += fob->chi[k] * de;
  if(x > 0.0){
    for(int k=0 ; k<fob->getNsize() ; k++) fob->chi[k] /= x;
  }

  /*** print average number of neutrons as function of TKE */
  FFPOutputnuTKE(fob->getKsize(),fob->nTKE);

  /*** print photon and neutron spectra */
  FFPOutputSpectrum(beohZeroCut(fob->spectrum[gammaray]),beohZeroCut(fob->chi),de,maxtemp,fob->spectrum,fob->chi);
}


/**********************************************************/
/*      Print Y(Z,A,TKE,dTKE)                             */
/**********************************************************/
void FFPOutputPairs(const double ycut, FissionFragmentPair *ffp)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  cout << "#   n   Zl   Al   Zh   Ah";
  cout << "    yield         TKE [MeV]     TXE [MeV]     dTXE [MeV]" << endl;

  double checksum = 0.0;
  int c = 0;
  for(int k=0 ; k<ffp->getN() ; k++){
    if(ffp->fragment[k].yield <= ycut) continue;

    cout << setw(5) << ++c;
    cout << setw(5) << ffp->fragment[k].getZl();
    cout << setw(5) << ffp->fragment[k].getAl();
    cout << setw(5) << ffp->fragment[k].getZh();
    cout << setw(5) << ffp->fragment[k].getAh();

    cout << setw(14) << ffp->fragment[k].yield;
    cout << setw(14) << ffp->fragment[k].ek;
    cout << setw(14) << ffp->fragment[k].ex;
    cout << setw(14) << ffp->fragment[k].sigma << endl;

    checksum +=  ffp->fragment[k].yield;
  }

  cout << "# Sum                    " << setw(14) << checksum << endl;
}


/**********************************************************/
/*      Print TXE Distribution                            */
/**********************************************************/
void FFPOutputTXEDistribution(FissionFragmentPair *ffp)
{
  FFPFindARange(ffp);

  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  for(int ah=ahmin ; ah<=ahmax ; ah++){

    double sum = 0.0;
    double txe = 0.0;
    for(int k=0 ; k<ffp->getN() ; k++){
      if((unsigned int)ah == ffp->fragment[k].getAh()){
        sum += ffp->fragment[k].yield;
        txe += ffp->fragment[k].yield * ffp->fragment[k].ex;
      }
    }
    if(sum > 0.0) txe /= sum;

    cout << setw(5) << ah;
    cout << setw(5) << ffp->af - ah;
    cout << setw(14) << txe << endl;
  }
}


/**********************************************************/
/*      Print TKE Distribution                            */
/**********************************************************/
void FFPOutputTKEDistribution(const int type, FissionFragmentPair *ffp)
{
  if(type <= 0 || type > 3) return;

  const int np = 300;
  double *px = new double [np];

  FFPFindARange(ffp);

  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  /* distribution for each Ah */
  for(int ah=ahmin ; ah<=ahmax ; ah++){

    double y = 0.0;
    if(type == 3){
      for(int k=0 ; k<ffp->getN() ; k++){
        if((unsigned int)ah == ffp->fragment[k].getAh()) y += ffp->fragment[k].yield;
      }
    }

    for(int i=0 ; i<np ; i++) px[i] = 0.0;

    /* for each Z,A */
    for(int k=0 ; k<ffp->getN() ; k++){
      if((unsigned int)ah != ffp->fragment[k].getAh()) continue;
      
      /* TKE distribution */
      double t0 = ffp->fragment[k].ek;
      double d0 = ffp->fragment[k].sigma;

      for(int i=0 ; i<np ; i++){
        double t1 = i + 0.5;
        px[i] += ffp->fragment[k].yield * exp(-(t0 - t1)*(t0 - t1) / (2.0*d0*d0));
      }
    }

    /* normalize TKE(A) distribution */
    double s = 0.0;
    for(int i=0 ; i<np ; i++) s += px[i];
    if(s > 0.0){
      for(int i=0 ; i<np ; i++) px[i] /= s;
    }

    /* average TKE(A) and variance */
    if(type == 1){
      double z = 0.0, v = 0.0;
      for(int i=0 ; i<np ; i++) z += (i + 0.5) * px[i];
      for(int i=0 ; i<np ; i++) v += (z - (i + 0.5)) * (z - (i + 0.5)) * px[i];

      cout << setw(5) << ah;
      cout << setw(14) << z;
      cout << setw(14) << sqrt(v) << endl;
    }

    /* histogram TKE distribution for a given A */
    else if(type == 2){
      for(int i=100 ; i<=250 ; i++){
        cout << setw(5) << ah;
        cout << setw(5) << i;
        cout << setw(14) << px[i] << endl;
      }
      cout << endl;
      cout << endl;
    }

    /* 3D data of Y(A,TKE) */
    else if(type == 3){
      for(int i=100 ; i<=250 ; i++){
        cout << setw(5) << ah;
        cout << setw(5) << i;
        cout << setw(14) << px[i] * y << endl;
      }
      cout << endl;
    }
  }

  delete [] px;
}


/**********************************************************/
/*      Print Gaussian Parameter and Chain Yield Y(A)     */
/**********************************************************/
void FFPOutputChainYield(const int n, const int ac, const int a0, double *s, double *d, double *f)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  cout << "#         Gauss1      Gauss2      Gauss3      GaussSym" << endl;
  cout << "# fract ";
  for(int k=0 ; k<4 ; k++) cout << setw(12) << f[k];
  cout << endl;
  cout << "# sigma ";
  for(int k=0 ; k<4 ; k++) cout << setw(12) << s[k];
  cout << endl;
  cout << "# delta ";
  for(int k=0 ; k<3 ; k++) cout << setw(12) << d[k];
  cout << endl;



  for(int i=0 ; i<n ; i++){
    double cy = 0.0;
    int a = i + a0;
    for(int k=0 ; k<4 ; k++){
      if(s[k] > 0.0){
        cy += gaussian(f[k],a - ac + d[k],s[k]);
        if(d[k] > 0.0) cy += gaussian(f[k],a - ac - d[k],s[k]);
      }
    }
    cout << setw(5) << a << setw(12) << cy << endl;
  }
  nl(); nl();
}


/**********************************************************/
/*      Print Individual CMS/LAB Spectrum                 */
/**********************************************************/
void FFPOutputIndividualSpectrum(ZAPair *zap, const int m, const double de, double *chiL, double *chiH, double *chi)
{
  const double emax = 30.0;

  cout << "#                     ";
  cout << setw(6) << zap->getZl();
  cout << setw(5) << zap->getAl();
  cout << setw(6) << zap->getZh();
  cout << setw(5) << zap->getAh();
  cout << "  Total    ";
  cout << " Cumulative" << endl;

  for(int k=0 ; k<=m ; k++){
    double e0 = (k > 0)  ? (k - 0.5)*de : 0;
    double e1 = (k + 0.5)*de;

    outVal(11,4,e0);
    outVal(11,4,e1);
    outVal(11,lowfilter(chiL[k]));
    outVal(11,lowfilter(chiH[k]));
    outVal(11,lowfilter(chiL[k] + chiH[k]));
    outVal(11,lowfilter(chi[k]));
    nl();

    if(e1 >= emax) break;
  }
  nl(); nl();
}


/**********************************************************/
/*      Print Average Multiplicity                        */
/**********************************************************/
void FFPOutputMultiplicity(double ytot, double *mult, double *eave, double *etot, double ep, double np)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  if(ytot == 0.0) return;

  ytot *= 0.5; // because ytot goes up to 2.0

  double x1 = etot[1]/ytot + ep;
  double x2 = mult[1]/ytot + np;
  double x3 = x1 / x2;

  outSectionHead("AVERAGE MULTIPLICITY AND ENERGY");
  if(np == 0.0){
    cout << "#                            gamma     neutron" << endl;
    cout << "# Total Energy[MeV]   ";
    cout << setw(12) << etot[0]/ytot << setw(12) << etot[1]/ytot << endl;
    cout << "# Average Energy[MeV] ";
    cout << setw(12) << eave[0]/ytot << setw(12) << eave[1]/ytot << endl;
    cout << "# Average Multiplicity";
    cout << setw(12) << mult[0]/ytot << setw(12) << mult[1]/ytot << endl;
  }
  else{
    cout << "#                            gamma     neutron +prefission" << endl;
    cout << "# Total Energy[MeV]   ";
    cout << setw(12) << etot[0]/ytot << setw(12) << etot[1]/ytot << setw(12) << x1 << endl;
    cout << "# Average Energy[MeV] ";
    cout << setw(12) << eave[0]/ytot << setw(12) << eave[1]/ytot << setw(12) << x3 << endl;
    cout << "# Average Multiplicity";
    cout << setw(12) << mult[0]/ytot << setw(12) << mult[1]/ytot << setw(12) << x2 << endl;
  }
  nl(); nl();
}


/**********************************************************/
/*      Print Multiplicity Distribution                   */
/**********************************************************/
void FFPOutputPnu(const int n, double *pnu)
{
  outSectionHead("NEUTRON MULTIPLICITY DISTRIBUTION");
  cout << "#  N  Pn(N)" << endl;

  double s = 0.0;
  for(int i=0 ; i<n ; i++){
    cout << setw(4) << i;
    outVal(12,pnu[i]); nl();

    s += i * pnu[i];
  }
  cout << "# av"; outVal(12,s); nl();

  nl(); nl();
}


/**********************************************************/
/*      Print Average Multiplicity                        */
/**********************************************************/
void FFPOutputYield(const int m, double *cy0, double *cy1, double **mult, double **eave)
{
  int a0 = 0, a1 = 0;
  for(int a=0 ; a<m ; a++){
    if((cy0[a] != 0.0) || (cy1[a] != 0.0)){ a0 = a-1 ; break; }
  }
  for(int a=m-1 ; a>=0 ; a--){
    if((cy0[a] != 0.0) || (cy1[a] != 0.0)){ a1 = a+1 ; break; }
  }

  outSectionHead("INDEPENDENT MASS CHAIN YIELD");
  cout << "#  A  preY(A)     postY(A)    pg(A)       pn(A)       <Eg>(A)     <En>(A)" << endl;
  for(int a=a0 ; a<=a1 ; a++){
    cout << setw(4) << a;
    outVal(12,cy0[a]);
    outVal(12,cy1[a]);
    outVal(12,mult[0][a]);
    outVal(12,mult[1][a]);
    outVal(12,eave[0][a]);
    outVal(12,eave[1][a]); nl();
  }
  nl(); nl();
}


/**********************************************************/
/*      Print nu(TKE                                      */
/**********************************************************/
void FFPOutputnuTKE(const int km, double *nTKE)
{
  int k0 = 0, k1 = 0;
  for(int k=1 ; k<km-1 ; k++){ if(nTKE[k] > 0.0){ k0 = k-1; break;} }
  for(int k=km-2 ; k>k0 ; k--){ if(nTKE[k] > 0.0){ k1 = k+1; break;} }


  outSectionHead("AVERAGE NUMBER OF NEUTRONS AS FUNCTION OF TKE");
  cout << "# TKE min   max[MeV]  nu" << endl;

  for(int k=k0 ; k<=k1 ; k++){
    outVal(10,2,(double)k);
    outVal(10,2,(double)k+1.0);
    outVal(12,nTKE[k]); nl();
  }
  nl(); nl();
}


/**********************************************************/
/*      Print Energy Spectrum                             */
/**********************************************************/
void FFPOutputSpectrum(const int m0, const int m1, const double de, const double t, double **spc, double *chi)
{
  int m = (m0 > m1) ? m0 : m1;

  outSectionHead("TOTAL PHOTON AND NEUTRON ENERGY SPECTRA");
  cout << "#      Emin       Emax"
       << "   Spectra/Fission/MeV"
       << "  Chi      ";
  if(t > 0.0) cout << " Maxw Ratio";
  cout << endl;
  cout << "#     [MeV]      [MeV]";

  cout << "      gamma    neutron     [/MeV]";
  if(t > 0.0){ cout <<"  T="; outVal(4,2,t); cout << "MeV"; }
  cout << endl;

  double x1[3], x2[3];
  for(int i=0 ; i<3 ; i++) x1[i] = x2[i] = 0.0;

  for(int k=0 ; k<=m ; k++){
    double e0 = (k > 0)  ? (k - 0.5)*de : 0;
    double e1 = (k + 0.5)*de;
    double e  = (e0 + e1) * 0.5;


    outVal(11,4,e0);
    outVal(11,4,e1);

    for(int p=0 ; p<2 ; p++) outVal(11,lowfilter(spc[p][k]));
    outVal(11,lowfilter(chi[k]));

    if(t > 0.0){
      double x = (k == 0) ? 0.25 * de : k*de;
      double p = 2/sqrt(PI*t*t*t) * exp(-x/t) * sqrt(x);
      outVal(11,lowfilter(chi[k] / p));
    }

    for(int p=0 ; p<2 ; p++){
      x1[p] += spc[p][k] * de;
      x2[p] += spc[p][k] * de * e;
    }
    x1[2] += chi[k] * de;
    x2[2] += chi[k] * de * e;
    nl();
  }
  nl();

  cout << "# Average Energy [MeV]";
  for(int i=0 ; i<3 ; i++){
    if(x1[i] > 0.0) outVal(11,x2[i]/x1[i]);
    else outVal(11,0.0);
  }

  nl(); nl();
}


/**********************************************************/
/*      Look for A-range in the Yield Data                */
/**********************************************************/
void FFPFindARange(FissionFragmentPair *ffp)
{
  /* find min/max Ah in FissionPair object */
  ahmin = 1000, ahmax = 0;
  for(int k=0 ; k<ffp->getN() ; k++){
    int ah = ffp->fragment[k].getAh();
    if(ah < ahmin) ahmin = ah;
    if(ah > ahmax) ahmax = ah;
  }
}


/**********************************************************/
/*      Gaussian Broadening Results by Given Resolution   */
/**********************************************************/
void FFPGaussianBroadening(const int m, const double dm, double *y0)
{
  double *y1 = new double [m];

  int d = (int)(3.0 * dm) + 1;
  for(int i0=0 ; i0<m ; i0++){

    double z = 0.0;
    y1[i0] = 0.0;
    for(int i1=i0-d ; i1<=i0+d ; i1++){
      if((i1 < 0) || (i1 >= m)) continue;
      double w = exp(-(i0 - i1)*(i0 - i1)/(dm*dm));
      z += w;
      y1[i0] += w * y0[i1];
    }
    y1[i0] = (z > 0.0) ? y1[i0]/z : 0.0;
  }
  for(int i=0 ; i<m ; i++) y0[i] = y1[i];

  delete [] y1;
}

