/******************************************************************************/
/*  dwcalc.cpp                                                                */
/*        DWBA calculation for weakly coupled states                          */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "structur.h"
#include "optical.h"
#include "omcalc.h"
#include "global.h"
#include "output.h"
#include "omoutput.h"
#include "etc.h"
#include "coupling.h"
#include "terminate.h"


static double  dwSpectroscopicAmplitude   (int, int);
static int     dwRadialIntegralZerorange  (int, int, int, int, double, int, CCdata *,
                                           complex<double> *, complex<double> *,
                                           complex<double> *, complex<double> *);
static void    dwMagneticSum              (int, int, int, int, int, int, int, int, int, int, int,
                                           double, complex<double>, complex<double> *);
static complex<double> dwOverlap          (int, double, complex<double> *, complex<double> *, complex<double> *);
static double  dwAngularDistribution      (int, int, int, double, int,
                                           double *, double *, CCdata *, complex<double> *);


extern double  *fact;

/**********************************************************/
/*     DWBA Calculation for Inelastic Scattering          */
/**********************************************************/
int dwbaCalc
(double        energy,        // CMS incident energy
 double        excitation,    // target excitation energy
 Pdata        *proj,          // particle data for projectile
 ZAnumber     *targ,          // target ZA
 double        mu,            // reduced mass for Targ+Proj system
 Direct       *dir,           // direct reaction data
 CrossSection *crx)           // calculated cross sections (output)
{
  Optical      omp[2];
  Potential    *pot;
  Wavefunc     *wfn = NULL;
  CCdata       cdt[2];
  LevelData    *lev = NULL;

  complex<double> *beta = NULL, *formfactor = NULL, *dwave0 = NULL, *dwave1 = NULL;

  int          l1max = 0,zzprod = 0;
  int          ndir0 = dir->ncc;

  if(ndir0 == 0) return(-1);
  if(energy <= 0.0) return(-1);
  if(dir->lev[ndir0].energy <= 0.0) return(-1);

  /*** no DWBA calculation for excited targets (to be improved) */
  if(excitation > 0) return(0);

//---------------------------------------
//     Memory Allocation

  try{
    pot = new Potential [3];
    for(int i=0 ; i<3 ; i++){
      pot[i].memalloc(MAX_POINTS);
    }
    wfn = new Wavefunc [2];
    for(int i=0 ; i<2 ; i++){
      wfn[i].memalloc(MAX_L,MAX_POINTS);
    }
    beta       = new complex<double> [MAX_LTRANS*MAX_L*9];
    dwave0     = new complex<double> [MAX_POINTS*MAX_L*3];
    dwave1     = new complex<double> [MAX_POINTS*MAX_L*3];

    lev = new LevelData [MAX_DIRECT];
  }
  catch(bad_alloc &e){
    message << "memory allocation error";
    cohTerminateCode("dwbaCalc");
  }

  zzprod = targ->getZ()*proj->za.getZ();

  lev[0].excitation = 0.0;
  lev[0].spin = abs(dir->lev[0].spin);

//---------------------------------------
//     Setup Energy and Optical Potential

  cdt[0].lev = &lev[0];
  omSetEnergy(energy, zzprod, mu, &lev[0]);

  omSetOmp(proj->omp,lev[0].energy,targ->getZ(),targ->getA(),proj->za.getZ(),proj->za.getA(),&omp[0]);

  omPotentialForm(zzprod, &omp[0], &lev[0], &pot[0]);

  /*** Store entrance channel distorted waves */
  dwStoreDistortedWave(proj->spin2,dir->nonloc,&cdt[0],&pot[0],&wfn[0],dwave0);

//---------------------------------------
//     For All Vibrational States

  int st = 0; // no spin flip
  int j0 = (int)abs((2.0*dir->lev[0].spin));
  int jc = 0;

  /*** Weak coupling model, core spin assumed, target spin - 1/2 */
  if((j0%2) != 0) jc = j0 - 1;

  int ndir = ndir0;
  for(int i=ndir0 ; i<MAX_DIRECT ; i++){

    if(dir->lev[i].energy <= 0.0) continue;

    /*** closed channel */
    if(dir->lev[i].energy >= energy) continue;

    /*** angular momentum transfer by difference in level spins */
    int lt = (int)fabs(abs(dir->lev[i].spin) - abs(dir->lev[0].spin));
    int jt = lt+st;

    lev[i].spin       = abs(dir->lev[i].spin);
    lev[i].excitation = dir->lev[i].energy;

    /*** Spectroscopic Amplitude */
    double spec_amp  = dwSpectroscopicAmplitude(st,proj->spin2)*dir->defdwba[i]/sqrt(2*lt+1.0);
    if(jc != 0) spec_amp *= (2*dir->lev[i].spin + 1.0)/( (j0+1.0)*(jc+1.0) );

    /*** Energy and Coulomb parameter for exit channel */
    cdt[1].lev = &lev[i];
    omSetEnergy(energy - dir->lev[i].energy, zzprod, mu, &lev[i]);

    /*** Exit channel optical potential geometry */
    omSetOmp(proj->omp,lev[i].energy,targ->getZ(),targ->getA(),proj->za.getZ(),proj->za.getA(),&omp[1]);

    omPotentialForm(zzprod, &omp[1], &lev[i], &pot[1]);
    pot[2].n_match = min(pot[0].n_match,pot[1].n_match);

    /*** Store exit channel distorted waves */
    dwStoreDistortedWave(proj->spin2,dir->nonloc,&cdt[1],&pot[1],&wfn[1],dwave1);

    /*** Collective form factor */
    int nint = dwCollectiveForm(spec_amp,&omp[0],&pot[2]);
    formfactor = pot[2].mean_field;

//---------------------------------------
//     Radial Integral and Scattering Amplitude

    if(nint){
      double coeff    = (j0+2*jt+1.0)/(j0+1.0) / (proj->spin2+1.0);

      l1max = dwRadialIntegralZerorange(lt,st,jt,nint,pot[0].width,proj->spin2,
                                        cdt,formfactor,beta,dwave0,dwave1);

      double tot = dwAngularDistribution(lt,l1max,MAX_ANGDIST,coeff,proj->spin2,crx->costh,
                                         crx->angdist[i],cdt,beta);
      crx->direct[i] = tot;
    }
    ndir++;
  }

//---------------------------------------
//     Output Cross Sections

  if(prn.xsection) outLevelExcite(ndir0,ndir,energy,excitation,lev,crx->direct);
  if(prn.angdist){
    /*** print out every ANGLE_STEP angle points */
    outAngularDistribution(2,ndir0,ndir,ANGLE_STEP,targ);
  }

//---------------------------------------
//     Free Allocated Memory

  delete [] beta;
  delete [] dwave0;
  delete [] dwave1;
  delete [] lev;

  return(l1max);
}


/**********************************************************/
/*     Distorted Wave Calculation                         */
/**********************************************************/
void dwStoreDistortedWave(const int spin2, double nonloc, CCdata *cdt,
                          Potential *pot, Wavefunc *wfn, complex<double> *dwave)
{
  complex<double>  *psave;
  complex<double> a,b,c,d;

  cdt->lmax = omExternalFunction(0,pot->rho_match,cdt->lev->coulomb,wfn);
  psave = wfn->internal;

  double spin = spin2/2.0;
  int i0   = spin2;
  int smax = i0+1;
  for(int l=0 ; l<=cdt->lmax ; l++){
    a = wfn->extderiv[l] / wfn->external[l] * pot->rho_match;
    b = conj(wfn->external[l])/ wfn->external[l];

    for(int s=0 ; s<smax ; s++){
      int ss = i0-2*s;
      if( (2*l+ss) < 0 ) continue;
      wfn->internal = &dwave[(3*l+s)*MAX_POINTS];

      double xj = l + (double)ss/2.0;
      if(xj < fabs(l-spin)) continue;
      omInternalFunction(pot->n_match,pot->width,cdt->lev->wavesq,(double)l,spin,xj,pot,wfn);
      d = omSmatrix(pot->n_match,pot->width,pot->rad_match,a,b,wfn);
      c = omNormalizationFactor(l,pot->n_match-3,cdt->lev->coulomb,cdt->lev->coulomb_scat0,d,wfn);

      omNormalization(pot->n_match,nonloc,c,pot,wfn);
    }
  }
  wfn->internal = psave;
}


/**********************************************************/
/*     Spectroscopic Amplitude                            */
/**********************************************************/
double dwSpectroscopicAmplitude(int st, int sa)
{
   double sfac=0.0;

   if(st == 0){
      sfac = sqrt(sa+1.0);
   }else if(st == 1){
      sfac = sqrt(3*(sa+1.0));
   }
   return(sfac);
}


/**********************************************************/
/*     Zero Range Radial Integral                         */
/**********************************************************/
int dwRadialIntegralZerorange(int lt, int st, int jt, int nint, double width,
                              int spin2, CCdata *cdt, complex<double> *form,
                              complex<double> *beta, complex<double> *dwave0, complex<double> *dwave1)
{
  complex<double> overlap, *wfn0,*wfn1; 
  int     i0,i1,j0,j1;
  double  x1,x2;

  int s0max = (i0=spin2)+1;
  int s1max = (i1=spin2)+1;
  int llt = lt*2;  
  st *= 2;
  jt *= 2;

  for(int k=0 ; k<MAX_LTRANS*MAX_L*9L ; k++) beta[k] = complex<double>(0.0,0.0);

  /*** Entrance channel wave function */
  int    l1maxt = 0;
  double y1 = sqrt((llt+1.0)*(st+1.0));

  for(int l0=0 ; l0<=cdt[0].lmax ; l0++){
    int ll0   = 2*l0;
    int l1min = abs(l0-lt);
    int l1max =     l0+lt ; if(l1max > cdt[1].lmax) l1max=cdt[1].lmax;
    double y2=sqrt(ll0+1.0)*y1;

    for(int s0=0 ; s0<s0max ; s0++){
      if((j0=ll0+(i0-2*s0)) < 0 ) continue;
      wfn0 = &dwave0[(3*l0+s0)*MAX_POINTS];

      /***  Exit channel wave function */
      for(int l1=l1min ; l1<=l1max ; l1+=2){
        int ll1 = 2*l1;
        if((x1=clebsh_gordan(ll1,llt,0,0,ll0)) == 0.0) continue;
        double z1 =(ll1+1) * ( ((l0-l1-lt)/2)%2==0 ? 1: -1);
        for(int s1=0 ; s1<s1max ; s1++){
          if( (j1=ll1+(i1-2*s1)) < 0 ) continue;
          if((x2=wigner_9j(jt,llt,st,j0,ll0,i0,j1,ll1,i1)) == 0.0) continue;
          double y3 =sqrt(j1+1.0)*y2;
          wfn1 = &dwave1[(3*l1+s1)*MAX_POINTS];

          overlap = dwOverlap(nint,width,form,wfn0,wfn1);
          dwMagneticSum(s0max,s1max,lt,jt,i0,i1,ll0,ll1,l1,j0,j1,z1*x1*x2*y3,overlap,beta);
        }
        if(l1 > l1maxt) l1maxt = l1;
      }
    }
  }

  return(l1maxt);
}


/**********************************************************/
/*     Magnetic momentum summation                        */
/**********************************************************/
void dwMagneticSum(int s0max,int s1max, int lt, int jt, int i0, int i1, 
                   int ll0, int ll1, int l1, int j0, int j1,
                   double z, complex<double> overlap, complex<double> *beta)
{
  double x1,x2,x3;
  int    m;

  for(int m0=0 ; m0<s0max ; m0++){
    int mm0 = i0-2*m0;
    if((x1=clebsh_gordan(ll0,i0,0,mm0,j0)) == 0.0) continue;
    for(int m1=0 ; m1<s1max ; m1++){
      int mm1 = i1-2*m1;
      for(int mp=0 ; mp<=lt ; mp++){
        int mm = 2*(m=mp+(mm1-mm0)/2);
        int m2 = abs(m);
        if(m2 > l1) continue;
        if((x2 = clebsh_gordan(ll1,i1,   -mm,       mm1,j1)) == 0.0) continue;
        if((x3 = clebsh_gordan( j1,jt,mm1-mm,mm-mm1+mm0,j0)) == 0.0) continue;
        double c = z*x1*x2*x3*sqrt((exp(fact[l1-m2]-fact[l1+m2])));
        unsigned long index = (unsigned int)MAX_L*(MAX_LTRANS*(3*m0+m1)+mp)+l1;
        beta[index] += overlap * c;
      }
    }
  }
}


/**********************************************************/
/*     Calculate Overlap Integral                         */
/**********************************************************/
complex<double> dwOverlap(int nint, double width, complex<double> *form, complex<double> *wfn0, complex<double> *wfn1)
{
  double f0 = width    /3.0;
  double f1 = width*2.0/3.0;
  double f2 = width*4.0/3.0;

  complex<double> x   = wfn0[nint] * wfn1[nint];
  complex<double> sum = x * form[nint] * f0;

  for(int i=1 ; i<nint ; i++){
    sum += wfn0[i] * wfn1[i] * form[i] * ((i%2==0) ? f1 : f2);
  }
  sum *= 2*SQRTPI;

  return(sum);
}


/**********************************************************/
/*     Angular distribution                               */
/**********************************************************/
double dwAngularDistribution(int lt, int l1max, int nth, double coeff, int spin2,
                             double *th, double *ad, CCdata *cdt, complex<double> *beta)
{
  complex<double> complex_Zero(0.0,0.0);
  int     i0,i1;

  int s0max = (i0=spin2)+1;
  int s1max = (i1=spin2)+1;

  double sum=0.0;
  for(int i=0 ; i<nth ; i++) ad[i]=0.0;
  for(int s0=0 ; s0<s0max ; s0++){
    int m0 =i0-2*s0;
    for(int s1=0 ; s1<s1max ; s1++){
      int m1 =i1-2*s1;
      for(int mp=0 ; mp<=lt ; mp++){
        int m  = mp+(m1-m0)/2;
        int m2 = abs(m);
        double x1 = (mp<0) ? ( (m2%2==0) ? 1 : -1 ) : 1;
        double x2 = (mp==0) ? 1 : 2;
        unsigned long index = MAX_L*(MAX_LTRANS*(3*s0+s1)+mp);

        for(int i=0 ; i<nth ; i++){
          complex<double> s(0.0,0.0);
          for(int l1=m2 ; l1<=l1max ; l1++){
            double plm = x1*assocLegendrePol(l1,m2,th[i]);
            s += beta[index+l1] * plm;
          }
          ad[i] += x2 * norm(s);
        }

        for(int l1=m2 ; l1<=l1max ; l1++){
          double x3 = 2.0/(2*l1+1.0)*( exp(fact[l1+m2]-fact[l1-m2]) );
          sum += x2*x3*norm(beta[index+l1]);
        }
      }
    }
  }
  double x1 = NORM_FACT*coeff/(PI4*PI4 * cdt[1].lev->energy * cdt[0].lev->energy)
             * cdt[1].lev->wave_number / cdt[0].lev->wave_number;

  for(int i=0 ; i<nth ; i++) ad[i] *= x1;
  sum *= x1*PI2;

  return(sum);
}

