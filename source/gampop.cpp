/******************************************************************************/
/*  gampop.cpp                                                                */
/*        gamma ray population                                                */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "structur.h"
#include "statmodel.h"
#include "levden.h"

static double  specTransisionGammaCont         (Statcalcmode, int, int, int, int, int,
                                                double, double, Nucleus *);
static double  specTransisionGammaDisc         (Statcalcmode, int, int, int, int, int,
                                                double, double, Nucleus *);

/**********************************************************/
/*      Gamma-ray Transition                              */
/**********************************************************/
double specTransitionGamma
( Statcalcmode  mode,              // calc. control for sum T/Moldauer/HF
  int           k0,                // index of parent excited state
  int           p0,                // parent state parity (-1/+1)
  int           j0,                // parent spate spin (doubled)
  double        **tg,              // transmission coeff. for gamma-ray
  double        q,                 // pop(parent)/Sigma T, if mode=HF
  Nucleus       *n,                // nucleus information
  double        *spc,              // gamma-ray spectrum
  double        *dlp)              // poplation increment
{
  /*** continuum to continuum */
  double x,y,sum = 0.0;
  int    kg;

  for(int k1=k0+1 ; k1<n->ncont ; k1++){

    x   = specTransisionGammaCont(mode,-1,k1,p0,j0,1,tg[E1][k1],q,n);
    x  += specTransisionGammaCont(mode, 1,k1,p0,j0,1,tg[M1][k1],q,n);
    x  += specTransisionGammaCont(mode, 1,k1,p0,j0,2,tg[E2][k1],q,n);
    x  += specTransisionGammaCont(mode,-1,k1,p0,j0,2,tg[M2][k1],q,n);
    x  += specTransisionGammaCont(mode,-1,k1,p0,j0,3,tg[E3][k1],q,n);

    if( (mode == hauser) || (mode == fluctuation) ){
      kg = k1-k0;
      y  = x/n->de;
      spc[kg] += y;
      dlp[k1] += y;
    }
    sum += x;
  }

  /*** continuum to level */
  for(int i=0 ; i<n->ndisc ; i++){

    if(n->lev[i].flag == 1) continue;

    double eg = n->excitation[k0] - n->lev[i].energy;
    if(eg < 0.0) continue;

    double ex = n->lev[i].energy;
    double a  = ldDensityParameter(ex,(double)n->za.getA(),&n->ldp);

    double te1 = gdrGammaTransmission(GL,E1,eg,a  ,ex );
    double tm1 = gdrGammaTransmission(SL,M1,eg,0.0,0.0);
    double te2 = gdrGammaTransmission(SL,E2,eg,0.0,0.0);
    double tm2 = gdrGammaTransmission(SL,M2,eg,0.0,0.0);
    double te3 = gdrGammaTransmission(SL,E3,eg,0.0,0.0);

    x   = specTransisionGammaDisc(mode,-1,i,p0,j0,1,te1,q,n);
    x  += specTransisionGammaDisc(mode, 1,i,p0,j0,1,tm1,q,n);
    x  += specTransisionGammaDisc(mode, 1,i,p0,j0,2,te2,q,n);
    x  += specTransisionGammaDisc(mode,-1,i,p0,j0,2,tm2,q,n);
    x  += specTransisionGammaDisc(mode,-1,i,p0,j0,3,te3,q,n);

    if( (mode == hauser) || (mode == fluctuation) ){
      kg = specFindEnergyBin(eg,n->de);
      if(kg >= 0){
        y  = x/n->de;
        spc[kg   ] += y;
        dlp[kg+k0] += y;
      }
    }
    sum += x;
  }

  return(sum);
}


/**********************************************************/
/*      Gamma-ray Transition : From Level to Levels       */
/**********************************************************/
double  specLevelTransisionGamma(Statcalcmode mode, int k0, int p0, int j0,
                                 double q, Nucleus *n)
{
  double x = 0.0;

  for(int k1=0 ; k1<k0-1 ; k1++){

    double eg = n->lev[k0].energy - n->lev[k1].energy;
    double ex = n->lev[k0].energy;
    double a  = ldDensityParameter(ex,(double)n->za.getA(),&n->ldp);

    double te1 = gdrGammaTransmission(GL,E1,eg,a  ,ex );
    double tm1 = gdrGammaTransmission(SL,M1,eg,0.0,0.0);
    double te2 = gdrGammaTransmission(SL,E2,eg,0.0,0.0);
    double tm2 = gdrGammaTransmission(SL,M2,eg,0.0,0.0);
    double te3 = gdrGammaTransmission(SL,E3,eg,0.0,0.0);

    x   = specTransisionGammaDisc(mode,-1,k1,p0,j0,1,te1,q,n);
    x  += specTransisionGammaDisc(mode, 1,k1,p0,j0,1,tm1,q,n);
    x  += specTransisionGammaDisc(mode, 1,k1,p0,j0,2,te2,q,n);
    x  += specTransisionGammaDisc(mode,-1,k1,p0,j0,2,tm2,q,n);
    x  += specTransisionGammaDisc(mode,-1,k1,p0,j0,3,te3,q,n);
  }

  return(x);
}

 
/**********************************************************/
/*      Gamma-ray Transition : Continuum to Continuum     */
/**********************************************************/
double specTransisionGammaCont(Statcalcmode mode, int pg, int k, int p0, int j0, int l,
                double tg, double q, Nucleus *n)
{
  double  x0=0.0,x1=0.0,dx=0.0,sum=0.0,wfc=1.0;
  int j1min = abs(j0-2*l);
  int j1max =     j0+2*l ;

  if(tg == 0.0) return(sum);

  if(j1max > 2*MAX_J-1) j1max=2*MAX_J-1;
  if(n->jmax < j1max/2) n->jmax=j1max/2;

  wfc = (mode == fluctuation) ? statMoldauer(-1,0,0,0,tg,0.0) : 1.0;

  for(int j1=j1min ; j1<=j1max ; j1+=2){
    if( (j0 == 0) && (j1 == 0) ) continue; // gamma decay prohibited for 0->0
    int idx = (j1-(int)(2.0*halfint(n->lev[0].spin)))/2;

    if(pg*p0 == 1){
      dx = n->density[k][idx].even*n->de;
      x0 = tg*dx;
      if(mode == wfactor){
        statWidthFluctuation(tg,dx,0.0);
      }else if( (mode == hauser) || (mode == fluctuation) ){
        n->pop[k][idx].even += (x1=x0*q*wfc);
      }
    }else{
      dx = n->density[k][idx].odd *n->de;
      x0 = tg*dx;
      if(mode == wfactor){
        statWidthFluctuation(tg,dx,0.0);
      }else if( (mode == hauser) || (mode == fluctuation) ){
        n->pop[k][idx].odd  += (x1=x0*q*wfc);
      }
    }
    if(mode == sumall) sum += x0;
    else if( (mode == hauser) || (mode == fluctuation) ) sum += x1;
  }

  return(sum);
}


/**********************************************************/
/*      Gamma-ray Transition : Continuum to Levels        */
/**********************************************************/
double specTransisionGammaDisc(Statcalcmode mode, int pg, int i, int p0, int j0, int l,
                double tg, double q, Nucleus *n)
{
  double x1=0.0,sum=0.0,wfc=1.0;
  double s1min = abs(j0/2.0-l);
  double s1max =     j0/2.0+l ;

  if(tg == 0.0) return(sum);
  if( (j0 == 0) && ((int)n->lev[i].spin == 0) ) return (sum);

  wfc = (mode == fluctuation) ? statMoldauer(-1,0,0,0,tg,0.0) : 1.0;

  if(n->lev[i].spin >=s1min && n->lev[i].spin<=s1max){
    if(pg == p0*n->lev[i].parity){
      if(mode == wfactor){
        statWidthFluctuation(tg,1.0,0.0);
      }else if( (mode == hauser) || (mode == fluctuation) ){
        n->lpop[i] += (x1=tg*q*wfc);
      }
      if(mode == sumall) sum += tg;
      else if( (mode == hauser) || (mode == fluctuation) ) sum += x1;
    }
  }

  return(sum);
}

