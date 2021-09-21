/******************************************************************************/
/*  photoabs.cpp                                                              */
/*        photon absorption cross section                                     */
/******************************************************************************/

#include <iostream>
#include <cmath>

using namespace std;

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "global.h"
#include "terminate.h"


static double photoGDR   (double, GDR *);
static double photoQD    (const double, const int, const int);

//  Fraction of QD, which goes into Pre-equilibrium calc
static double qdratio = 0.0;
double photoQDratio()
{ return qdratio; }

/**********************************************************/
/*     Gamma-ray induced entrance channel calculation     */
/*         GDR model and quasi deuteron break-Up          */
/**********************************************************/
double photoAbsorption(double eg, double kg, int targlev, Nucleus *n, Transmission *tin, GDR *gdr)
{
  if(eg <= 0.0){
    message << "incident gamma-ray energy " << eg << " negative or zero";
    cohTerminateCode("photoAbsorption");
  }

  for(int j=0 ; j<3*MAX_J ; j++) tin->tran[j] = 0.0;

  double te1 = gdrGammaTransmission(SL,E1,eg,n->ldp.a,n->excitation[0]);
  double tm1 = gdrGammaTransmission(SL,M1,eg,0.0,0.0);
  double te2 = gdrGammaTransmission(SL,E2,eg,0.0,0.0);

  /*** GDR absorption */
  double absgdr = 0.0;
  for(int i=0 ; i<MAX_GDR ; i++){
    absgdr += photoGDR(eg,&gdr[i]);
  }

  /*** recalculate SigR from transmissions */
  absgdr =(  te1 * 3
           + tm1 * 3
           + te2 * 5) * PI/kg/kg * NORM_FACT / 2.0;

  /*** quasi-deutron absorption */
  double absqd  = photoQD(eg,n->za.getZ(),n->za.getA());

  /*** add QD to GDR, and renormalize, is it OK ? */
  qdratio = (absgdr > 0.0) ? absqd/absgdr : 0.0;
  double f = 1.0 + qdratio;
  qdratio = absqd / (absqd + absgdr);  if(qdratio > 1.0) qdratio = 1.0;

  /*** insert transmission coefficients at L=1 and L=2 */
  tin->tran[3] = (te1+tm1)*f/2.0;
  tin->tran[6] = te2*f;
  tin->lmax    = 2;

  /*** set max J = I + 2 */
  n->jmax = n->lev[targlev].spin + 2.0;

  return(absgdr*f);
}


/**********************************************************/
/*      GDR Absorption                                    */
/**********************************************************/
double photoGDR(double eg, GDR *gdr)
{
  double e = gdr->getEnergy();
  double w = gdr->getWidth();
  double s = gdr->getSigma();

  double a = eg * eg * w * w;
  double b = e * e - eg * eg;

  return(s * a / (b * b + a));
}


/**********************************************************/
/*      Quasi Deuteron Absorption                         */
/*      Chadwick et al., PRC44, 814 (1991)                */
/**********************************************************/
double photoQD(const double eg, const int z, const int a)
{
  const double deuteron_bind = 2.22452;
  const double levinger_const = 6.5;
  double sigd = 0.0, pbf = 0.0;

  double e2 = eg * eg;
  double e3 = e2 * eg;
  double e4 = e3 * eg;

  if(eg > deuteron_bind) sigd = 61.2*( pow((eg-2.224),1.5) )/e3;

  if(eg < 20.0)       pbf = exp(-73.3   /eg);
  else if(eg > 140.0) pbf = exp(-24.2348/eg);
  else                pbf = 0.083714 - 0.0098343*eg + 4.1222e-04*e2 - 3.4762e-06*e3 + 9.3537e-09*e4;

  double x = (double)(a-z)*(double)z/(double)a;

  return(levinger_const * x * sigd * pbf);
}

