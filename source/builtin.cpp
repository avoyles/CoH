/******************************************************************************/
/*  builtin.cpp                                                               */
/*        built in global optical model potentials                            */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>

using namespace std;

#include "omplib/omplib.h"
#include "etc.h"

static int          findname    (string *, char *, int);
static unsigned int add_library (int, double, int, int, int, int, Optical *);

static const int N_OMPLIB = 23;
static string lib_name[N_OMPLIB] =
  {"Wilmore", "Becchetti", "Rapaport", "Walter", "CH89", 
   "Koning", "Soukhovitskii", "Kunieda",
   "Perey", "Schwandt", "Madland",
   "Lemos", "Nolte", "Avrigeanu", "Avrigeanu2009", "Avrigeanu2014", "TALYS-A",
   "Bojowald", "An", "Han",
   "bound", "test", "spline"};

static const int ADD_CUSTOM_OMP = 17;
static string add_name[ADD_CUSTOM_OMP] =
  {"scratch",
   "SmithA120","SmithA136","SmithA415",
   "Flap22","Young-Am","Young-Pu","Young-Re","modSoukhovitskii","Soukhovitskii2005",
   "Dave1p","WLH1","WLH2","WLH3","WLH4","WLH5"};


/***********************************************************/
/*      Search for Given OMP Names in the Library          */
/***********************************************************/
unsigned int find_omp(string str)
{
   unsigned int k = findname(lib_name,&str[0],N_OMPLIB);
   if(k==0){
        if( (k  = findname(add_name,&str[0],ADD_CUSTOM_OMP))!=0 ) k +=N_OMPLIB;
   }
   if(k==0) cerr << "ERROR     : OMP name [ " << str << " ] not found" << endl;
   return(k);
}


/***********************************************************/
/*      Compare OMP Names Provided With Library            */
/***********************************************************/
int findname(string *list, char *name, int n)
{
  int  i;

  for(i=0 ; i<n ; i++){
    if(!strcmp(name,list[i].c_str())) break;
  }
  if(i==n)  return(0);
  else      return(i+1);
}


/***********************************************************/
/* Potential form bit field code description               */
/*   bit: 76543210  :  0 Imag  Volume  (WS)                */
/*                     1       Surface (deriv. WS)         */
/*                     2       Surface (gaussian)          */
/*                     3       undef.                      */
/*                     4 Real  Volume  (WS)                */
/*                     5       Surface (deriv. WS)         */
/*                     6       Surface (gaussian)          */
/*                     7       undef.                      */
/*   example                                               */
/*        00010001 = 0x0011 : Real WS + Imag WS            */
/*        00010010 = 0x0012 : Real WS + Imag dWS           */
/*        00010011 = 0x0013 : Real WS + Imag WS +dWS       */
/***********************************************************/
/***********************************************************/
/*      Select Optical Potential Parameter                 */
/***********************************************************/
unsigned int omp_library(int id, int zt, int at, int zi, int ai, 
                         double e, Optical *omp)
{
  unsigned int potfm;
  switch(id){
  case  1: potfm = OMPWilmoreHodgson(      at,   ai,zi,omp); break;
  case  2: potfm = OMPBecchettiGreenlees(  at,zt,ai,zi,omp); break;
  case  3: potfm = OMPRapaport(          e,at,zt,ai,zi,omp); break;
  case  4: potfm = OMPWalterGuss(        e,at,zt,ai,zi,omp); break;
  case  5: potfm = OMPChapelHill89(      e,at,zt,ai,zi,omp); break;
  case  6: potfm = OMPKoningDelaroche(   e,at,zt,ai,zi,omp); break;
  case  7: potfm = OMPSoukhovitskii(     e,at,zt,   zi,omp); break;
  case  8: potfm = OMPKunieda(           e,at,zt,ai,zi,omp); break;
  case  9: potfm = OMPPerey(               at,zt,ai,zi,omp); break;
  case 10: potfm = OMPMadlandSchwandt(   e,at,zt,ai,zi,omp); break;
  case 11: potfm = OMPMadlandYoung(      e,at,zt,      omp); break;
  case 12: potfm = OMPLemos(                     ai,zi,omp); break;
  case 13: potfm = OMPNolte(               at,zt,ai,zi,omp); break;
  case 14: potfm = OMPAvrigeanu(         e,at,zt,ai,zi,omp); break;
  case 15: potfm = OMPAvrigeanu2009(     e,at,zt,ai,zi,omp); break;
  case 16: potfm = OMPAvrigeanu2014(     e,at,zt,ai,zi,omp); break;
  case 17: potfm = OMPTALYS_alpha(       e,at,zt,ai,zi,omp); break;
  case 18: potfm = OMPBojowald(          e,at,zt,ai,zi,omp); break;
  case 19: potfm = OMPAnHaixia(            at,zt,ai,zi,omp); break;
  case 20: potfm = OMPHanYinlu(            at,zt,ai,zi,omp); break;
  case 21: potfm = OMPBoundState(                   zi,omp); break;
  case 22: potfm = OMPtest(                         zi,omp); break;
  case 23: potfm = OMPspline(            e,at,zt,ai,zi,omp); break;
  default: potfm = 0                                       ; break;
  }

  if(id>N_OMPLIB) potfm = add_library(id-N_OMPLIB,e,zt,at,zi,ai,omp);

  return(potfm);
}


unsigned int add_library(int id, double e, int zt, int at, int zi, int ai, Optical *omp)
{
  unsigned int potfm;
  switch(id){
  case  1: potfm = OMPscratch(           e,at,zt,ai,zi,omp); break;
  case  2: potfm = OMPSmithA120(         e,at,zt,      omp); break;
  case  3: potfm = OMPSmithA136(         e,at,zt,      omp); break;
  case  4: potfm = OMPSmithA415(           at,zt,      omp); break;
  case  5: potfm = OMPFlap22(            e,at,zt,      omp); break;
  case  6: potfm = OMPYoung_Am(          e,at,zt,      omp); break;
  case  7: potfm = OMPYoung_Pu(          e,at,zt,      omp); break;
  case  8: potfm = OMPYoung_Re(          e,at,zt,      omp); break;
  case  9: potfm = OMPModSoukhovitskii(  e,at,zt,   zi,omp); break;
  case 10: potfm = OMPSoukhovitskii2005( e,at,zt,   zi,omp); break;
  case 11: potfm = OMPDave1p(              at,zt,      omp); break;
//  case 12: potfm = OMPWLH1(              e,at,zt,ai,zi,omp); break;
//  case 13: potfm = OMPWLH2(              e,at,zt,ai,zi,omp); break;
//  case 14: potfm = OMPWLH3(              e,at,zt,ai,zi,omp); break;
//  case 15: potfm = OMPWLH4(              e,at,zt,ai,zi,omp); break;
//  case 16: potfm = OMPWLH5(              e,at,zt,ai,zi,omp); break;
  default: potfm = 0; break;
  }

  return(potfm);
}


/*************************************************/
/*    Dispersion Contribution to Surface Term    */
/*    J.M. Quesada et al. PRC 67, 067601 (2003)  */
/*************************************************/
double OMPdeltaSurface(double ex, double e0, double as, double bs, double cs)
{
  const int m = 2;

  double ep = ex + e0;
  double em = ex - e0;

  double resp = - ep*ep/(ep*ep + bs*bs);
  double resm =   em*em/(em*em + bs*bs);

  complex<double> vs(0.0,0.0), p, z;
  for(int j=1 ; j<=m ; j++){
    p = (j==1) ? complex<double>(0.0,bs) : complex<double>(0.0,-bs);
    z = ex/m * p*(2.0*p + ep - em)/((p + e0)*(p + ep)*(p - em));
    complex<double> z1 = -p*cs;
    vs += z * exp(z1) * expintE1(z1);
  }
  vs = as*(vs - resp * exp( cs*ep) * expint(-cs*ep)
              - resm * exp(-cs*em) * expint( cs*em))/PI;

  return(vs.real());
}


/*************************************************/
/*    Dispersion Contribution to Volume Term     */
/*************************************************/
double OMPdeltaVolume(double ex, double e0, double av, double bv)
{
  const int m = 2;

  double ep = ex + e0;
  double em = ex - e0;

  double resp = - ep*ep/(ep*ep + bv*bv);
  double resm =   em*em/(em*em + bv*bv);

  complex<double> vv(0.0,0.0), p, z;
  for(int j=1 ; j<=m ; j++){
    p = (j==1) ? complex<double>(0.0,bv) : complex<double>(0.0,-bv);
    z = ex/m * p*(2.0*p + ep - em)/((p + e0)*(p + ep)*(p - em));
    vv += z * log(-p);
  }
  vv = -av*(vv + resp * log(abs(ep)) + resm * log(abs(em)))/PI;
  
  return(vv.real());
}


