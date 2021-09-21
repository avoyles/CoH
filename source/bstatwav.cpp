/******************************************************************************/
/*  bstatwav.cpp                                                              */
/*        bound state wave functions                                          */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "structur.h"
#include "optical.h"
#include "dsd.h"

static void        bstateWsWave (double, double, double, Bstate *, Potential *);
static inline bool is_magic     (int);

/**********************************************************/
/*      Bound Wave Calculation for Shell Model            */
/**********************************************************/
int bstateWavefunc(Pdata *proj, ZAnumber *targ, double eb, double mu,
                   double nonloc, double bmax, double fshift, Optical *omp, Potential *pot, Bstate *bst)
{
  int       maxorbit;
  double    u,v;

  int m = (proj->pid == proton) ? targ->getZ() : targ->getN();

  /*** Quantum numbers of all bound states */
  int nfermi = nilssonFillOrbit(proj->pid,m);

  for(int i=0 ; i<MAX_SPLEVEL ; i++){
    if( (nilssonLevel(proj->pid,i,nfermi,targ->getA(),&bst[i],eb,bmax,fshift)) < 0 ){
      bst[i].n = bst[i].l = bst[i].j2 = 0;  bst[i].bind = 0.0;
    }
  }

  /*** Look for highest level having negative binding energy */
  maxorbit = MAX_SPLEVEL;
  for(int i=0 ; i<MAX_SPLEVEL ; i++){
    if(bst[i].bind>0.0){
      maxorbit = i;
      break;
    }
  }

  /*** BCS calculation, find Fermi surface */
  double delta   = bcsEnergyGap(proj->pid,targ->getZ(),targ->getA());
  double efermi  = (bst[nfermi+1].bind + bst[nfermi].bind)/2.0;
  efermi = bcsGapEquation(m,efermi,delta,bst);

  for(int i=0 ; i<MAX_SPLEVEL ; i++){
    if(is_magic(m)){
      u = (i < nfermi) ? 0.0 : 1.0;
    }else{
      v = bcsVfactor(efermi,bst[i].bind,delta);
      u = sqrt(1-v*v);
    }
    bst[i].w = u;
  }


  /*** solve Schroedinger equation for each bound state */
  for(int i=0 ; i<maxorbit ; i++) bstateWsWave(mu,omp->volume.real(),nonloc,&bst[i],pot);
/*
  for(int i=0 ; i<MAX_SPLEVEL ; i++){
    if(bst[i].j2==0) break;
    printf("%3d %3d %3d %3d % 7.3f % 10.3f",i,
           bst[i].n, bst[i].l,bst[i].j2, bst[i].bind, bst[i].w*bst[i].w);
    if(i==nfermi) cout << "  <= Nf";
    cout << endl;
  }
*/
  return(pot->n_match);
}


/**********************************************************/
/*      Bound Wave Calculation (Woods Saxon form)         */
/*      for fixed (n,l,j), adjust V                       */
/**********************************************************/
void bstateWsWave(double mu, double v0, double nonloc, Bstate *bst,Potential *pot)
{
  if(bst->bind >= 0.0){
    for(int i=0 ; i<=pot->n_match ; i++) bst->wave[i]=0.0;
    return;
  }
  double k2    = bst->bind*mu*2.0*AMUNIT/VLIGHTSQ/HBARSQ;
  double alpha = mu*nonloc*nonloc*AMUNIT/VLIGHTSQ/HBARSQ/2.0;

  v0 = bstateFindDepth(bst->n,bst->l,bst->j2, k2,v0,alpha,bst->wave,pot);

// for(i=0;i<=pot->n_match;i++) printf("%d %e\n",i,bst->wave[i]); printf("\n\n");
}



/**********************************************************/
/*     Check if Magic Number                              */
/**********************************************************/
bool is_magic(int m)
{
  const int magic_number[8] = {2,8,20,28,50,82,126,186};
  bool magic = false;
  for(int k=0 ; k<8 ; k++){
    if(m==magic_number[k]){
      magic = true;
      break;
    }
  }
  return(magic);
}
