/******************************************************************************/
/*  popinit.cpp                                                               */
/*        initial population in the first compound                            */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"

#undef DEBUG_POPCHECK

/**********************************************************/
/*      Store initial population at the given excitation  */
/**********************************************************/
void statStoreInitPopulation(int spin0, int targlev, int targid, double coef, Transmission *tin)
{
  int ig    = (int)(2.0*halfint(ncl[0].lev[0].spin));   // half-int ofset for CN
  int it    = (int)(2.0*ncl[targid].lev[targlev].spin); // target state spin
  int pt    = ncl[targid].lev[targlev].parity;          // target state parity
  double g  = coef /( (it+1.0)*(spin0+1.0) );           // spin factor without (2J+1)

  /*** Loop over CN J and Parity*/
  for(int j0=ig ; j0<=ncl[0].jmax*2+ig ; j0+=2){
    int   jdx = (j0-ig)/2;
    for(int p0=-1 ; p0<=1 ; p0+=2){

      /*** For entrance channel Tlj coupled to JP*/
      for(int lp0=0 ; lp0<=tin->lmax*2 ; lp0+=2){
        if(parity(lp0) != p0*pt) continue;

        for(int sp0=spin0 ; sp0>=-spin0 ; sp0-=2){
          int jp0 = lp0 + sp0;
          if(jp0<0) continue;

          if( abs(jp0-it)>j0 || j0>(jp0+it) ) continue;

          double tlj  = tin->tran[tj_index(lp0,sp0,spin0)];
          if(tlj==0.0) continue;

          if(p0==1)  ncl[0].pop[0][jdx].even += (j0+1.0)*g*tlj;
          else       ncl[0].pop[0][jdx].odd  += (j0+1.0)*g*tlj;
        }
      }
    }
  }

#ifdef DEBUG_POPCHECK
  cout.setf(ios::fixed, ios::floatfield);
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  double s = 0.0;
  for(int j=0 ; j<20 ; j++) {
    s += ncl[0].pop[0][j].even + ncl[0].pop[0][j].odd;
    cout << setw(10) << j;
    cout << setw(12) << ncl[0].pop[0][j].even;
    cout << setw(12) << ncl[0].pop[0][j].odd;
    cout << endl;
  }
  cout << setw(12) << s << endl;
#endif
}


/**********************************************************/
/*      Adjust Initial Population to Match Reaction CS    */
/**********************************************************/
void statAdjustInitPopulation(double sigr)
{
  if(sigr == 0.0) return;

  double f = 0.0;
  for(int j=0 ; j<=ncl[0].jmax ; j++) f += ncl[0].pop[0][j].even + ncl[0].pop[0][j].odd;
  if(f == 0.0) return;

  f = sigr/f;
  for(int j=0 ; j<=ncl[0].jmax ; j++){
    ncl[0].pop[0][j].even *= f;
    ncl[0].pop[0][j].odd  *= f;
  }
  f=0.0;
  for(int j=0;j<=ncl[0].jmax;j++)  f += ncl[0].pop[0][j].even + ncl[0].pop[0][j].odd;
} 
