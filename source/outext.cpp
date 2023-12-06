/******************************************************************************/
/*  outext.cpp                                                                */
/*        non-standard output                                                 */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "outext.h"
#include "nucleus.h"

static const bool gammaray_with_conversion = false;

static inline void extHead (const char *section){
  cout << "#" << section << endl;
}

static inline void extHead (const char *section, const double x){
  cout << "#" << section << setprecision(4) << setw(12) << x << endl;
}

static inline void extHead (const char *section, ZAnumber za){
  cout << "#" << section;
  cout << setw(3) << za.getZ() << setw(4) << za.getA() << endl;
}


/**********************************************************/
/*      Cumulative Levels and Spin Distributions          */
/**********************************************************/
void extCumulativeLevels(Nucleus *n)
{
  int kmax = 0;
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    if(ncl[p].ndisc>kmax) kmax = ncl[p].ndisc;
  }
  if(kmax == 0) return;

  double *sdist[MAX_CHANNEL];
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    sdist[j] = new double [MAX_J];
    for(int l=0 ; l<MAX_J ; l++) sdist[j][l] = 0.0;
  }

  /*** heading for discrete levels */
  extHead("levels");
  cout << "#         ";
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    cout << "      " << setw(3) << ncl[p].za.getZ() << setw(5) << ncl[p].za.getA();
  }
  cout << endl;

  /*** for all discrete states */
  for(int k=0 ; k<kmax ; k++){
    cout << setw(11) << k+1;
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(!n->cdt[j].status) continue;
      int p = n->cdt[j].next;
      if( (k!=0 && ncl[p].lev[k].energy==0.0) || (ncl[p].ndisc<=k) ){
        cout << "  ?       ?   ";
      }else{
        char s = (ncl[p].lev[k].parity < 0) ? '-' : '+';
        cout.setf(ios::fixed, ios::floatfield);
        cout << setw(8) << setprecision(4) << ncl[p].lev[k].energy;
        cout << setw(5) << setprecision(1) << ncl[p].lev[k].spin << s;

        int l = (int)ncl[p].lev[k].spin;
        sdist[j][l] += 1.0;
      }
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;


  /*** normalize spin distributions */
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    double s = 0.0;
    for(int l=0 ; l<MAX_J ; l++) s += sdist[j][l];
    if(s > 0.0){
      for(int l=0 ; l<MAX_J ; l++) sdist[j][l] /= s;
    }
  }

  int lmax = 0;
  for(int l=MAX_J-1 ; l>=0 ; l--){
    bool nonzero = false;
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(sdist[j][l] > 0.0) nonzero = true;
    }
    if(nonzero){ lmax = l+1; break; }
  }

  /*** heading for spin distributions */
  extHead("spindist");
  cout << "#         J";
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    cout << "      " << setw(3) << ncl[p].za.getZ() << setw(5) << ncl[p].za.getA();
  }
  cout << endl;

  for(int l=0 ; l<=lmax ; l++){
    cout << setw(11) << l;
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(!n->cdt[j].status) continue;
      cout << setw(14) << setprecision(6) << sdist[j][l];
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;


  for(int j=0 ; j<MAX_CHANNEL ; j++) delete [] sdist[j];
}


/**********************************************************/
/*      Level Densities in Continuum                      */
/**********************************************************/
void extContinuumDensity(Nucleus *n)
{
  int kmax = 0;

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    if(ncl[p].ntotal > kmax) kmax = ncl[p].ntotal;
  }
  if(kmax == 0) return;

  for(int j=0 ; j<MAX_CHANNEL ; j++){

    if(!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    if(ncl[p].ntotal <= 1) continue;

    /*** heading */
    extHead("density", ncl[p].za);

    /*** for all continuum bins */
    cout << "# Ex[MeV]  ";
    cout << " Total[/MeV] ";

    double j0 = (((int)(ncl[p].lev[0].spin*2.0 + 0.001) % 2) == 1) ? 0.5 : 0.0;

    cout.setf(ios::fixed, ios::floatfield);
    for(int j=0 ; j<=ncl[p].jmax ; j++){
      cout << setprecision(1) << setw(13) << j + j0;
    }
    cout << endl;

    cout.setf(ios::scientific, ios::floatfield);
    for(int k=0 ; k<ncl[p].ntotal ; k++){

      if(ncl[p].excitation[k] <= 0.0) continue;

      double sum = 0.0;
      for(int j=0 ; j<=ncl[p].jmax ; j++){
        sum +=  ncl[p].density[k][j].even + ncl[p].density[k][j].odd;
      }

      cout << " " << setprecision(3) << setw(10) << ncl[p].excitation[k];
      cout << setprecision(6) << setw(13) << sum;
      for(int j=0 ; j<=ncl[p].jmax ; j++){
        cout << setprecision(6) << setw(13) << ncl[p].density[k][j].even + ncl[p].density[k][j].odd;
      }
      cout << endl;
    }
    cout << endl;
  }
}


/**********************************************************/
/*      Fission Level Densities                           */
/**********************************************************/
void extFissionDensity(const int nc)
{
  extHead("fissiondensity");
  cout.setf(ios::scientific, ios::floatfield);

  for(int i=0 ; i<nc ; i++){
    if(!ncl[i].fissile) continue;

    /*** heading */
    cout << "#";
    cout << setw(3) << ncl[i].za.getZ() << setw(4) << ncl[i].za.getA() << endl;

    cout << " Energy    Barrier 1    Barrier 2    Barrier 3 " << endl;

    /*** for all continuum bins */
    for(int k=0 ; k<ncl[i].ntotal ; k++){
      cout << setprecision(3) << setw(10) << ncl[i].excitation[k];

      for(int h=0 ; h<MAX_HUMP ; h++){
        double sum = 0.0;
        for(int j=0 ; j<MAX_J ; j++) sum +=  ncl[i].fission->barrier[h].density[k][j].even
                                           + ncl[i].fission->barrier[h].density[k][j].odd;
        cout << setprecision(6) << setw(13) << sum;
      }
      cout << endl;
    }
    cout << endl;
    cout << endl;
  }

  cout << endl;
}


/**********************************************************/
/*      Discrete Level and Gamma-Ray Branching Ratio      */
/**********************************************************/
void extGammaBranch()
{
  int    fs[MAX_GAMMA_BRANCH], f;
  double br[MAX_GAMMA_BRANCH], b;
  double gr[MAX_GAMMA_BRANCH], g;

  if(gammaray_with_conversion) extHead("gammaraywithconversion");
  else extHead("gammaray");

  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(6);


  for(int c=0 ; c<MAX_CHANNEL ; c++){

    if(!ncl[0].cdt[c].status) continue;
    int p = ncl[0].cdt[c].next;

    if(ncl[p].ndisc == 0) continue;

    /*** ID for CN and number of discrete levels included */
    cout << setw(5) << c << setw(4) << ncl[p].ndisc;
    cout << setw(3) << ncl[p].za.getZ() << setw(4) << ncl[p].za.getA() << endl;

    for(int i=0 ; i<ncl[p].ndisc ; i++){
      cout << setw(9) << i << setw(4) << ncl[p].lev[i].ngamma
           << setw(13) << ncl[p].lev[i].energy;

      /*** sort the gamma branches by decreasing order */
      for(int m=0 ; m<ncl[p].lev[i].ngamma ; m++){
        fs[m] = ncl[p].lev[i].fstate[m];
        br[m] = ncl[p].lev[i].branch[m];
        gr[m] = ncl[p].lev[i].gratio[m];
      }

      for(int j=0 ; j<ncl[p].lev[i].ngamma ; j++){
        int k = j;
        for(int m=j ; m<ncl[p].lev[i].ngamma ; m++){
          if(fs[m] > fs[k]) k=m;
        }
        f = fs[j];  fs[j] = fs[k];  fs[k] = f;
        b = br[j];  br[j] = br[k];  br[k] = b;
        g = gr[j];  gr[j] = gr[k];  gr[k] = g;
      }

      /*** print final state and branching ratio */
      for(int m=0 ; m<ncl[p].lev[i].ngamma ; m++){
        if(gammaray_with_conversion){
          cout << setw(4) << fs[m] << setw(13)<< br[m] << setw(13)<< gr[m];
        }
        else{
          cout << setw(4) << fs[m] << setw(13)<< br[m];
        }
      }
      cout << endl;
    }
  }
  cout << endl;
}


/**********************************************************/
/*      All Discrete Gamma Lines                          */
/**********************************************************/
void extGammaLine(const double elab, double *spc, Nucleus *n, const double gw)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(6);

  int kmax = (int)(n->max_energy/n->de) + 1; // max number of spectrum bin
  int gmax = kmax + 3 * (int)(gw/n->de);     // when Gauss width is given, add 3 sigma

  double *gcnt = new double [gmax+1];

  for(int k=0 ; k<=kmax ; k++){
    gcnt[k] = spc[k];

    double x0 = (k-0.5)*n->de; if(x0 < 0.0) x0 = 0.0;
    double x1 = (k+0.5)*n->de;

    /*** subtract discrete gamma lines */
    for(int i=0 ; i<gml.getN() ; i++){
      double eg = abs(gml.line[i].energy);
      if((x0 <= eg) && (eg < x1)) gcnt[k] -= gml.line[i].production/n->de;
    }
    /*** if negative or too small, set zero */
    if(gcnt[k] < 1e-10) gcnt[k] = 0.0;
  }

  extHead("photonproduction",elab);

  /*** not broadening */
  if(gw == 0.0){
    cout << "# discrete  " << setw(11) << gml.getN() << endl;
    cout << "#       Z   A Energy       Production" << endl;
    for(int i=0 ; i<gml.getN() ; i++){
      cout << setw(5) << i+1;
      cout << setw(4) << gml.line[i].za.getZ();
      cout << setw(4) << gml.line[i].za.getA();
      cout << setprecision(6);
      cout << setw(13) << abs(gml.line[i].energy);
      cout << setw(13) << gml.line[i].production << endl;
    }

    cout << endl;
    cout << endl;
    cout << "# continuum " << setw(11) << kmax+1 << setprecision(4) << setw(11) << n->de << endl;
    cout << "#       Emin         Emax       Total        Continuum" << endl;
    for(int k=0 ; k<=kmax ; k++){
      double x0 = (k-0.5)*n->de; if(x0 < 0.0) x0 = 0.0;
      double x1 = (k+0.5)*n->de;

      cout << setw(5) << k+1;
      cout << setprecision(4) << setw(13) << x0 << setw(13) << x1;
      cout << setprecision(6) << setw(13) << spc[k] << setw(13) << gcnt[k] << endl;
    }
  }
  /*** when broadening width given */
  else{
    /*** first, boaden the continuum part */
    specGaussianBroadening(gmax, n->de, gw, gcnt);

    /*** add discrete level contributions */
    for(int i=0 ; i<gml.getN() ; i++){
      specGaussianBroadening(gmax, n->de, gw, gcnt, abs(gml.line[i].energy), gml.line[i].production);
    }

    cout << "# continuum " << setw(11) << gmax+1 << setprecision(4) << setw(11) << n->de << endl;
    cout << "#       Emin         Emax       Total        Continuum" << endl;
    for(int k=0 ; k<=gmax ; k++){
      double x0 = (k-0.5)*n->de; if(x0 < 0.0) x0 = 0.0;
      double x1 = (k+0.5)*n->de;

      cout << setw(5) << k+1;
      cout << setprecision(4) << setw(13) << x0 << setw(13) << x1;
      cout << setprecision(6) << setw(13) << spc[k] << setw(13) << gcnt[k] << endl;
    }
  }

  delete [] gcnt;
}


/**********************************************************/
/*      Channel Decay Widths                              */
/**********************************************************/
void extDecayWidth(const int j, double *w)
{
  extHead("decaywidth");

  cout.setf(ios::fixed, ios::floatfield);
  cout << setw(5) << setprecision(1) << j/2.0;

  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(5);
  for(int i=0 ; i<MAX_CHANNEL+1 ; i++) cout << setw(13) << w[i];
  cout << endl;
}


/**********************************************************/
/*      Gamma Strength Function                           */
/**********************************************************/
void extGammaStrengthFunction(const int n, const double de, double **tg)
{
  extHead("gstrengthfunc");

  cout << "#   Eg[MeV]       E1[mb]       M1[mb]       E2[mb]       M2[mb]       E3[mb]" << endl;

  double c0 = HBARSQ * VLIGHTSQ * PI * PI * NORM_FACT / PI2; // convert transmission to photo-absorption

  for(int k=1 ; k<n ; k++){
    double eg = k * de;
    double e2 = eg * eg;
    double c1 = c0 * 3.0 / e2;
    double c2 = c0 * 5.0 / e2;
    double c3 = c0 * 7.0 / e2;

    cout.setf(ios::fixed, ios::floatfield);
    cout << setprecision(4) << setw(11) << eg;

    cout.setf(ios::scientific, ios::floatfield);
    cout << setprecision(5);
    cout << setw(13) <<  tg[k][E1] * c1;
    cout << setw(13) <<  tg[k][M1] * c1;
    cout << setw(13) <<  tg[k][E2] * c2;
    cout << setw(13) <<  tg[k][M2] * c2;
    cout << setw(13) <<  tg[k][E3] * c3;
    cout << endl;
  }
}

