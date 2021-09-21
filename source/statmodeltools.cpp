/******************************************************************************/
/*  statmodeltools.cpp                                                        */
/*        some utilitiy functions to calcualte statistical model              */
/******************************************************************************/

#include <iostream>
#include <cmath>

using namespace std;

#include "structur.h"
#include "nucleus.h"
#include "statmodel.h"


/**********************************************************/
/*      Find a Bin Index for Given Excitation             */
/**********************************************************/
int    specFindEnergyBin(double e, double de)
{
  double eps = 1.0e-5;

  if(fabs(e) <= eps) return(0);
  if( (e < 0.0) || (e > (MAX_ENERGY_BIN-1)*de) ) return(-1);

  if(e < de*0.5) return(0);

  /*** Add first energy bin, which has a half width,
       so that zero-energy emission is mapped onto k=0 */ 
  e += 0.5*de;

  /*** To avoid round-off error, a small number is added */
  int k = (int)((e+eps)/de);
  if( (de*k <= (e+eps)) && (e < de*(k+1)) ) return (k);

  cout << " out of range " << k << " " << de*k <<  " " << e << " " << de*(k+1) << endl;
  return(-1);
}


/**********************************************************/
/*      Sum Up All Particle Emission Spectra              */
/**********************************************************/
void specCumulativeSpectra(const int csize, const int nsize, double **spc, Nucleus *n)
{
  int cs = crx.getNspectra();
  int ns = crx.getNspecbin();

  for(int id=0 ; id<csize ; id++){
    if(n->cdt[id].status){
      for(int k=0 ; k<nsize ; k++){
        if((id < cs) && (k < ns)) crx.spectra[id][k] += spc[id][k];
      }
    }
  }
}


/**********************************************************/
/*      Gaussian Broadening Results by Given Resolution   */
/**********************************************************/
void specGaussianBroadening(const int m, const double de, const double dm, double *y0)
{
  double *y1 = new double [m];

  int d = 3 * (int)(dm / de) + 1;

  for(int i0=0 ; i0<m ; i0++){

    double z = 0.0;
    y1[i0] = 0.0;
    for(int i1=i0-d ; i1<=i0+d ; i1++){
      if((i1 < 0) || (i1 >= m)) continue;
      double e0 = i0 * de;
      double e1 = i1 * de;
      double w = exp(-(e0-e1)*(e0-e1)/(2*dm*dm));
      z += w;
      y1[i0] += w * y0[i1];
    }
    y1[i0] = (z > 0.0) ? y1[i0]/z : 0.0;
  }

  for(int i=0 ; i<m ; i++) y0[i] = y1[i];

  delete [] y1;
}


/**********************************************************/
/*      Add Gaussian Broadened Line to Continuum          */
/**********************************************************/
void specGaussianBroadening(const int m, const double de, const double dm, double *y0, const double x1, const double y1)
{
  int d = 3 * (int)(dm / de) + 1;
  double *y2 = new double [2*d+1];

  int i0 = -1;
  for(int i1=0 ; i1<m-1 ; i1++){
    double e0 = i1 * de;
    double e1 = e0 + de;
    if((e0 <= x1) && (x1 < e1)){
      i0 = i1;
      break;
    }
  }
  if(i0 < 0) return;

  double z = 0.0;
  for(int k = 0 ; k<=2*d ; k++){
    y2[k] = 0.0;

    int i1 = i0 - d + k;
    if((i1 < 0) || (i1 >= m)) continue;

    double e1 = i1 * de;
    double w = exp(-(x1-e1)*(x1-e1)/(2*dm*dm));
    z += w;
    y2[k] = w * y1;
  }

  if(z > 0.0){
    z *= de;
    for(int k = 0 ; k<=2*d ; k++) y2[k] = y2[k] / z;
  }

  for(int k = 0 ; k<=2*d ; k++) y0[i0-d+k] += y2[k];

  delete [] y2;
}

