/******************************************************************************/
/*  gdr.cpp                                                                   */
/*        global GDR parameter library                                        */
/******************************************************************************/

#include <iostream>
#include <cmath>

using namespace std;

#include "gdr.h"


/**********************************************************/
/*      E1                                                */
/*      M. Herman et al. EMPIRE-3.2 Malta Manual          */
/*      Eq. (4.108) -- (4.116)                            */
/**********************************************************/
void gdrE1(double mass, double beta, GDR *g)
{
  double energy = (49.336+7.34*beta)*pow(mass,-0.2409);
  double width  = 0.3 * energy;
  double sigma  = 10.6 * mass / width;
  string XL     = "E1";
  g->setGDR(XL,energy,width,sigma);
}

void gdrE1DoubleHump0(double mass, double beta, GDR *g)
{
  if(beta <= 0.064){
    gdrE1(mass,beta,g);
  }
  else{
    double energy = 50.0*pow(mass,-0.232) * exp(-0.946*beta);
    double width  = (0.282-0.263*beta) * energy;
    double sigma  = 3.48  * mass / width;
    string XL     = "E1";
    g->setGDR(XL,energy,width,sigma);
  }
}

void gdrE1DoubleHump1(double mass, double beta, GDR *g)
{
  if(beta <= 0.064){
    g->clear();
  }
  else{
    double energy = 50.0*pow(mass,-0.232);
    double width  = (0.35 -0.14 *beta) * energy;
    //double sigma  = 8.26  * mass / width;
    double sigma  = 1.464 * pow(mass,4.0/3.0) / width;
    string XL     = "E1";
    g->setGDR(XL,energy,width,sigma);
  }
}


/**********************************************************/
/*      M1                                                */
/**********************************************************/
void gdrM1(double mass, GDR *g)
{
  double energy = 41.0/pow(mass,1/3.0);
  double width  = 4.0;
  double sigma  = 1.0;
  string XL     = "M1";
  g->setGDR(XL,energy,width,sigma);
}


/**********************************************************/
/*      M1 Scissors                                       */
/**********************************************************/
void gdrM1scissors(double mass, double beta, GDR *g)
{
  double energy = 80.0*abs(beta)/pow(mass,1/3.0);
  double width  = 1.5;
  double sigma  = 42.4*beta*beta / width;
  string XL     = "M1";
  g->setGDR(XL,energy,width,sigma);
}


/**********************************************************/
/*      E2                                                */
/**********************************************************/
void gdrE2(double z, double mass, GDR *g)
{
  double a3 = 1./pow(mass,1/3.0);
  double energy = 63.0*a3;
  double width  = 6.11 - 0.012*mass;
  double sigma  = 1.5e-4 * z * z * energy * energy*a3 / width;
  string XL     = "E2";
  g->setGDR(XL,energy,width,sigma);
}


/**********************************************************/
/*      M2                                                */
/**********************************************************/
void gdrM2(GDR *gm1, GDR *g)
{
  const double ratio = 0.0008;

  double energy = gm1->getEnergy();
  double width  = gm1->getWidth();
  double sigma  = gm1->getSigma() * ratio;
  string XL     = "M2";
  g->setGDR(XL,energy,width,sigma);
}


/**********************************************************/
/*      E3                                                */
/**********************************************************/
void gdrE3(GDR *ge2, GDR *g)
{
  const double ratio = 0.0008;

  double energy = ge2->getEnergy();
  double width  = ge2->getWidth();
  double sigma  = ge2->getSigma() * ratio;
  string XL     = "E3";
  g->setGDR(XL,energy,width,sigma);
}

