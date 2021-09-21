/******************************************************************************/
/*  ecloutput.cpp                                                             */
/*        print out calculated exclusive spectra and DDX                      */
/******************************************************************************/

#include <iostream>
#include <iomanip>

using namespace std;

#include "structur.h"
#include "nucleus.h"
#include "eclipse.h"
#include "output.h"

//static inline void   eclHead (const char *section){
//  cout << "#" << section << endl;
//}


/**********************************************************/
/*      Top Part                                          */
/**********************************************************/
void eclOutHead(const int n, const int c, double elab)
{
  cout << setprecision(6) << setiosflags(ios::scientific);
  cout << "# ECLIPSE    "
       << setw(13) << elab
       << setw(13) << n
       << setw(13) << c    << endl;
}


/**********************************************************/
/*      Ouput Individual Spectrum                         */
/**********************************************************/
void eclOutSpectra(const int nm, const int km, const int cm, double de, double **np, EXSpectra *dat)
{
  cout << setprecision(6) << setiosflags(ios::scientific);

  for(int n=0 ; n<nm ; n++){

    /*** find non-zero element from the high-side */
    int  km2  = eclFindKmax(km,cm,dat[n].spec);

    if(km2 <= 0) continue;

    cout << "# ";
    outZA(&ncl[n].za);

    for(int c=0 ; c<cm ; c++) cout << setw(13) << np[n][c];
    cout << endl;

    for(int k=0 ; k<km2 ; k++){
      cout << setw(13) << k*de;
      for(int c=0 ; c<cm ; c++) cout << setw(13) << dat[n].spec[c][k];
      cout << setw(13) << dat[n].glin[k];
      if(n < cm)  cout << setw(13) << dat[n].sbin[k];
      cout << endl;
    }
    cout << endl;
    cout << endl;
    cout << endl;
  }
}


/**********************************************************/
/*      Ouput Total Spectra                               */
/**********************************************************/
void eclTotalSpectra(const int nm, const int km, const int cm, double de, EXSpectra *dat)
{
  double *sum;
  sum = new double [cm];

  /*** sum all the partial spectra to the first section in dat */
  for(int n=1 ; n<nm ; n++){
    for(int k=0 ; k<km ; k++){
      for(int c=0 ; c<cm ; c++) dat[0].spec[c][k] += dat[n].spec[c][k];
    }
  }

  /*** find non-zero element from the high-side */
  int  km2  = eclFindKmax(km,cm,dat[0].spec);

  if(km2 > 0){
    for(int c=0 ; c<cm ; c++){
      sum[c] = 0.0;
      for(int k=0 ; k<km2 ; k++) sum[c] += dat[0].spec[c][k];
    }

    cout << "# Total Spectra" << endl;

    for(int k=0 ; k<km2 ; k++){
      cout << setw(13) << k*de;
      for(int c=0 ; c<cm ; c++) cout << setw(13) << dat[0].spec[c][k];
      cout << endl;
    }
    cout << "# Sum        ";
    for(int c=0 ; c<cm ; c++) cout << setw(13) << sum[c];
    cout << endl;
    cout << endl;
    cout << endl;
  }

  delete [] sum;
}


/**********************************************************/
/*      Write Legendre Coefficient, Header Part           */
/**********************************************************/
void eclOutNucleusHead(const int n, const int cm, double **np)
{
  static char p[] = {'g','n','p','a','d','t','h'};
  int j = 0;

  cout << "# Nucleus    ";
  cout <<" ";
  for(int i=1 ; i<7 ; i++){
    j = (i >= cm) ? 0 : (int)np[n][i];
    cout << setw(1) << p[i] << setw(1)<< j;
  }
  cout << setw(13) << n << endl;
}


/**********************************************************/
/*      Write Legendre Coefficient, Sub-Header Part       */
/**********************************************************/
void eclOutChannelHead(const int c, double mp)
{
  int nsec = 0;

  if(c == 0) nsec = (mp > 0.0) ? 2 : 0;
  else       nsec = (mp > 0.0) ? 1 : 0;

  cout << "# Particle   ";
  cout << setw(13) << c
       << setw(13) << nsec
       << setw(13) << mp << endl;
}


/**********************************************************/
/*      Write Legendre Coefficient                        */
/**********************************************************/
void eclOutLegCoeff(const int km, const int nleg, double *ep, double **cleg)
{
  if(km <= 0){
    cout << "# Continuum  " << setw(13) << 0 << setw(13) << 0 << endl;
  }
  else{
    int k0 = 1, k1 = km-2, k2 = km;
    /*** cut low side zeros, keep the first two points */
    for(int k=2 ; k<km ; k++){
      bool p = false;
      for(int l=0 ; l<nleg ; l++) if(cleg[l][k] != 0.0) p = true;
      if(p){
        k0 = k-1;
        break;
      }
    }

    /*** high side check */
    for(int k=km-3 ; k>=k0 ; k--){
      bool p = false;
      for(int l=0 ; l<nleg ; l++) if(cleg[l][k] != 0.0) p = true;
      if(p){
        k1 = k+1;
        break;
      }
    }
    k2 = k1 - k0 + 3;

    cout << "# Continuum  " << setw(13) << k2 << setw(13) << nleg << endl;

    cout << setw(13) << ep[0];
    for(int l=0 ; l<nleg ; l++) cout << setw(13) << cleg[l][0];
    cout << endl;

    for(int k=k0 ; k<=k1 ; k++){
      cout << setw(13) << ep[k];
      for(int l=0 ; l<nleg ; l++) cout << setw(13) << cleg[l][k];
      cout << endl;
    }

    cout << setw(13) << ep[km-1];
    for(int l=0 ; l<nleg ; l++) cout << setw(13) << cleg[l][km-1];
    cout << endl;
  }
  cout << endl;
  cout << endl;
}


/**********************************************************/
/*      Write Discrete Gamma Lines                        */
/**********************************************************/
void eclOutGammaLine(const int ng, const int nleg, double *e, double *x)
{
  cout << "# Discrete   " << setw(13) << ng << setw(13) << nleg << endl;
  for(int j=0 ; j<ng ; j++){
    cout << setw(13) << e[j] << setw(13) << x[j]<< endl;
  }
  cout << endl;
  cout << endl;
}


/**********************************************************/
/*      Print Out Decay Probability Table                 */
/**********************************************************/
void eclOutPtable(const int n, const int m, int *ktot, int **cidx, double ****ptbl)
{
  cout.setf(ios::scientific, ios::floatfield);

  for(int c0=0 ; c0<n ; c0++){
    for(int k0=0 ; k0<=ktot[c0] ; k0++){

      cout << setw(5) << c0 << setw(4) << k0 << setw(4) << ktot[c0];
      for(int j=0 ; j<m ; j++) cout << setw(13) <<  cidx[c0][j];
      cout << endl;

      for(int k1=k0 ; k1<=ktot[c0] ; k1++){
        cout << setw(13) << k1;
        for(int c1=0 ; c1<m ; c1++){
          cout << setw(13) << ptbl[c0][k0][c1][k1];
        }
        cout << endl;
      }
    }
  }
}


/**********************************************************/
/*      Look For Non-Zero Component                       */
/**********************************************************/
int eclFindKmax(const int km, const int cm, double **dat)
{
  bool kend = false;
  int  km2  = km-1;

  for(int k=km-1 ; k>0 ; k--){
    for(int c=0 ; c<cm ; c++){
      if(dat[c][k] > 0.0){
        kend = true;
        break;
      }
    }
    if(kend){
      km2 = k+2;
      break;
    }
  }
  return( (kend) ? km2 : -1 );
}
