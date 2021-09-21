/******************************************************************************/
/*  stattrans.cpp                                                             */
/*        set-up transmission coefficients                                    */
/******************************************************************************/
#include <iostream>
#include <cmath>

using namespace std;

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"
#include "omcalc.h"
#include "levden.h"
#include "masstable.h"
#include "parameter.h"


static const double ecms_cutoff = 1e-9; // lowest emission energy 1 meV

//#define ModifiedLorentzian
//#define StandardLorentzian

#undef DEBUG_TRANS
#ifdef DEBUG_TRANS
#include <cstdio>
#endif


/**********************************************************/
/*      Store Particle Continuum Transmission into Array  */
/**********************************************************/
void statStoreContinuumTransmission(int ip, double ex, Pdata *p, Transmission **tc)
{
  CrossSection cx;

  /*** Clear transmission array */
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[j][k].memclear();
  }

  /*** For all particle decay channels */
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    if (!ncl[ip].cdt[j].status) continue;

    double sp = ncl[ip].cdt[j].binding_energy;
    int    id = ncl[ip].cdt[j].next;
    double mu = ncl[id].mass * p[j].mass / (ncl[id].mass + p[j].mass);
    double f  = parmGetFactor(parmTJ,ncl[ip].za,(Particle)j);

    for(int k=0 ; k<ncl[id].ntotal ; k++){
      tc[j][k].ecms = ex - ncl[id].excitation[k] - sp;

      /*** For the lowest transition, represented by 1/4 mesh point */
      if(k==0){
        tc[j][k].ecms += (ncl[id].excitation[0] - ncl[id].excitation[1])/4.0;
      }
      if(tc[j][k].ecms < ecms_cutoff) continue;
      
      tc[j][k].lmax = omCalc(tc[j][k].ecms, &p[j], &ncl[id].za, mu, tc[j][k].tran, &cx);
      tc[j][k].sigr = cx.reaction * f;

      for(int l=0 ; l<=3*tc[j][k].lmax ; l++) tc[j][k].tran[l] *= f;
#ifdef DEBUG_TRANS
      printf("CONT %3d  %3d%3d %5d % 10.5f % 10.5f % 11.4e % 11.4e % 11.4e\n",
             j,k,id,tc[j][k].lmax,tc[j][k].ecms,ncl[id].excitation[k],
             tc[j][k].tran[0],tc[j][k].tran[3],tc[j][k].tran[4]);
#endif
    }

    /*** Account for quater-size mesh for k=0 */
    for(int l=0 ; l<=3*tc[j][0].lmax ; l++) tc[j][0].tran[l] *= 0.25;
  }
}


/**********************************************************/
/*      Store Particle Discrete Transmission into Array   */
/**********************************************************/
void statStoreDiscreteTransmission(int ip, double ex, Pdata *p, Transmission **td)
{
  CrossSection cx;

  /*** Clear transmission array */
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    for(int k=0 ; k<MAX_LEVELS ; k++) td[j][k].memclear();
  }

  /*** For all particle decay channels */
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    if (!ncl[ip].cdt[j].status) continue;

    double sp = ncl[ip].cdt[j].binding_energy;
    int    id = ncl[ip].cdt[j].next;
    double mu = ncl[id].mass * p[j].mass / (ncl[id].mass + p[j].mass);
    double f  = parmGetFactor(parmTJ,ncl[ip].za,(Particle)j);

    /*** Transition to discrete levels */
    for(int k=0 ; k<ncl[id].ndisc ; k++){
      td[j][k].ecms = ex - ncl[id].lev[k].energy - sp;
      if(td[j][k].ecms < ecms_cutoff) continue;

      td[j][k].lmax = omCalc(td[j][k].ecms, &p[j], &ncl[id].za, mu, td[j][k].tran, &cx);
      td[j][k].sigr = cx.reaction * f;

      for(int l=0 ; l<=3*td[j][k].lmax ; l++) td[j][k].tran[l] *= f;
#ifdef DEBUG_TRANS
      printf("DISC %3d  %3d%3d %5d % 10.5f % 10.5f % 11.4e % 11.4e % 11.4e % 11.4e % 11.4e % 11.4e % 11.4e\n",
             j,k,id,td[j][k].lmax,td[j][k].ecms,ncl[id].lev[k].energy,
             td[j][k].tran[0],
             td[j][k].tran[3],td[j][k].tran[4],
             td[j][k].tran[6],td[j][k].tran[7],
             td[j][k].tran[9],td[j][k].tran[10]);
#endif
    }
  }
}


/**********************************************************/
/*  Store Gamma-ray Transmissions into 2-dim Array        */
/**********************************************************/
void statStoreGammaTransmission(int k0, double **tg, Nucleus *n)
{
  /*** Clear transmission array */
  for(int i=0 ; i<MAX_MULTIPOL ; i++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tg[i][k] = 0.0;
  }

  /*** Neutron separation energy for temperature dependent Gamma */
  double sn = 0.0;
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if( n->cdt[id].pid == neutron ) sn = n->cdt[id].binding_energy;
  }
  if(sn == 0.0){
    bool checkmass = true;
    double mx = mass_excess(n->za.getZ(),n->za.getA()-1,&checkmass);
    sn = (checkmass) ? mx - n->mass_excess + ENEUTRON : 0.0;
  }

  for(int k1=k0+1 ; k1<n->ncont ; k1++){
    double eg = n->excitation[k0] - n->excitation[k1];
    double ex = (eg <= sn) ? sn - eg : 0.0;
    double a  = ldDensityParameter(ex,(double)n->za.getA(),&n->ldp);

#ifdef ModifiedLorentzian
    tg[E1][k1] = gdrGammaTransmission(ML,E1,eg,a  ,ex );
#elif defined StandardLorentzian
    tg[E1][k1] = gdrGammaTransmission(SL,E1,eg,a  ,ex );
#else
    tg[E1][k1] = gdrGammaTransmission(GL,E1,eg,a  ,ex );
#endif
    tg[M1][k1] = gdrGammaTransmission(SL,M1,eg,0.0,0.0);
    tg[E2][k1] = gdrGammaTransmission(SL,E2,eg,0.0,0.0);
    tg[M2][k1] = gdrGammaTransmission(SL,M2,eg,0.0,0.0);
    tg[E3][k1] = gdrGammaTransmission(SL,E3,eg,0.0,0.0);

#ifdef DEBUG_TRANS
    if(k0 == 0){
      printf(" %3d %10.6f %10.6f %10.6f %11.4e %11.4e %11.4e\n",
           k1,eg,ex, ((ex < 0.0) ? 0.0 : sqrt(ex/a)),
           tg[E1][k1]/(PI2*pow(eg,3.0)),
           tg[M1][k1]/(PI2*pow(eg,3.0)),
           tg[E2][k1]/(PI2*pow(eg,5.0)));
    }
#endif
  }
}
