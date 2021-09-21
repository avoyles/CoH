/******************************************************************************/
/*  ffpyield.cpp                                                              */
/*        Calculate Y(Z,A,TKE) distribution                                   */
/******************************************************************************/

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "masstable.h"
#include "wahl.h"
#include "ffp.h"
#include "terminate.h"

static const unsigned char DEBUG_MASS = 0x01;
static const unsigned char DEBUG_LIST = 0x02;
static const unsigned char DEBUG_NORM = 0x04;
static const unsigned char DEBUG_CHRG = 0x08;
static const unsigned char DEBUG_ENRG = 0x10;
static const unsigned char DEBUG_ATKE = 0x20;

//static unsigned char DEBUG_ID = DEBUG_MASS | DEBUG_LIST;
static unsigned char DEBUG_ID = 0x0;

static void FFPDebug (unsigned char, FissionFragmentPair *);

/* data for fission yield */
static double *tke, *dtke;                    // TKE(A) and its width
static double *chain_yield, **charge_dist;    // Y(A), Y(Z,A)
static int    *charge_first;                  // smallest Z number for Y(Z,A) data

static void   FFPZAList (FissionFragmentPair *);
static void   FFPYield (FissionFragmentPair *);
static void   FFPTXEDistribution (FissionFragmentPair *);
static void   FFPAllocateMemory (FissionFragmentPair *);
static void   FFPDeleteAllocated (FissionFragmentPair *);


/**********************************************************/
/*      Mass-Chain Yield by Gaussians                     */
/**********************************************************/
void FFPMassChainYieldParameters(const bool sf, string data, double *sigma, double *delta, double *fract, FissionFragmentPair *ffp)
{
  /*** internal Gaussian mass distribution based on experimental data */
  if(data == "internal"){
    FFPGaussianParameters(sf,ffp->zf,ffp->af,ffp->ex,ffp->bn,sigma,delta,fract);
  }
  /*** mass distribution by Wahl, set nu= 0 */
  else if(data == "Wahl"){
    WahlChainYield(ffp->zf,ffp->af,ffp->ex,sigma,delta,fract);
  }

  /*** re-normalize fractions, just in case */
  double s = 2.0*(fract[0] + fract[1] + fract[2]) + fract[3];
  double d = 2.0 - s;
  for(int i=0 ; i<4 ; i++) fract[i] += d * fract[i]/s;

  /*** copy parameters to FissionFragmentProduct object */
  ffp->resetN();
  for(int i=0 ; i<4 ; i++){
    ffp->setGauss(sigma[i], delta[i],fract[i]);
    if(delta[i] != 0.0) ffp->setGauss(sigma[i],-delta[i],fract[i]);
  }
}


/**********************************************************/
/*            Fragment Yield from Raw Data File           */
/**********************************************************/
void FFPReadRawYield(string filename, FissionFragmentPair *ffp)
{
  /*** Read the raw primary fragment yield data from formatted ASCII file */

  unsigned int maxA = 0;   // The maximum fragment mass
  unsigned int minA = 400; // The minimum fragment mass
  double tke_file[minA];   // Read in TKE(A)
  double dtke_file[minA];  // Read in dTKE(A)

  /*** allocate data arrays */
  FFPAllocateMemory(ffp);

  // File exists so continue
  ifstream infile(filename.c_str(),ifstream::in);
  if(!infile){
    message << "File " << filename << " does not exist!";
    cohTerminateCode("FFPReadRawYield");
  }

  double am = ffp->af * 0.5;

  // Read raw data
  string str;
  int k = 0; // nuclide counter
  while(getline(infile,str)){
    if(str[0] == '#') continue;

    unsigned int zl;         // The light fragment charge
    unsigned int al;         // The light fragment mass
    double y;                // The fragment pair yield
    double ke;               // The fragment pair total kinetic energy
    double dke;              // The fragment pair distribution in kinetic energy

    // Read data from file
    istringstream ss(str);
    ss >> zl >> al >> y >> ke >> dke;

    // Automatically skip the heavy fragment if it is provided
    if((double)al > am) continue;

    // Calculate split partner
    unsigned int zh = ffp->zf - zl;
    unsigned int ah = ffp->af - al;

    // Calculate Q value for partition
    double mcn = mass_excess(ffp->zf,ffp->af) + ffp->ex;
    double q   = mcn - (mass_excess(zl,al) + mass_excess(zh,ah));

    // Consistency check for total excitation energy
    if(q - ke > 0.0){
      // Save fragment information
      ffp->fragment[k].setPair(zl,al,zh,ah);
      ffp->fragment[k].qval  = q;      // Set the Q-value
      ffp->fragment[k].yield = y;      // Set the fragment yield
      ffp->fragment[k].ex = q - ke;    // Set the TXE
      ffp->fragment[k].ek = ke;        // Set the TKE
      tke_file[al] = tke_file[ah] = ke;
      dtke_file[al] = dtke_file[ah] = dke;
      if (al <= minA) { minA = al; }
      if (ah >= maxA) { maxA = ah; }
      k++;
    }
  }
  infile.close();

  // Set the total number of nuclides (assuming blank line in file there should be one less)
  ffp->setN(k-1);
  ffp->setZARange(minA, 20);

  // Set the number of Gaussians to zero
  ffp->setG(0);

  // Set TKE(A) from the read in data file and some generic systematics for width
  for(int i=0 ; i<ffp->nmass ; i++){
    int a = i + ffp->mass_first; // Offset which gives fragment mass, a
    tke[i] = tke_file[a];
    dtke[i] = dtke_file[a];    
  }
  
  // Finally, calculate width of TXE, assuming it is proportional to TKE
  for(int k=0 ; k<ffp->getN() ; k++){
    // Average TKE
    double tke0 = tke[ffp->fragment[k].getAl() - ffp->mass_first];
    // Set the width of TXE (this comes from error propagation of the <TKE>(A) - see Shin Okumura's paper)
    ffp->fragment[k].sigma = ffp->fragment[k].ek * dtke[ffp->fragment[k].getAl() - ffp->mass_first] / tke0;
  }

  // Renormalize the yields
  double s = 0.0;
  for(int k=0 ; k<ffp->getN() ; k++) s += ffp->fragment[k].yield;
  for(int k=0 ; k<ffp->getN() ; k++) ffp->fragment[k].yield /= s;
  
  // For debugging purposes print out yield averaged TKE so we know we're on the right track before the long calculation starts
  cout << "# <TKE> = " << FFPAverageTKE(ffp) << endl;

  // Debugging 
  // FFPDebug(DEBUG_LIST, ffp);
  // FFPDebug(DEBUG_ENRG | DEBUG_CHRG | DEBUG_ATKE, ffp);

  FFPDeleteAllocated(ffp);
}


/**********************************************************/
/*      Fission Fragment Pair and Yield Calculation       */
/**********************************************************/
void FFPGeneratePairs(const bool sf, const double tke_input, FissionFragmentPair *ffp, double *fzn)
{
  /*** allocate data arrays */
  FFPAllocateMemory(ffp);

  /*** mass chain yield by Gaussians */
  for(int i=0 ; i<ffp->nmass ; i++){
    chain_yield[i] = 0.0;
    int a = i + ffp->mass_first;
    for(int k=0 ; k<ffp->getNg() ; k++){
      if(ffp->getGaussWidth(k) > 0.0){
        double da = a - ffp->ac + ffp->getGaussCenter(k);
        chain_yield[i] += gaussian(ffp->getGaussFraction(k),da,ffp->getGaussWidth(k));
      }
    }
  }
  if(DEBUG_ID & DEBUG_MASS) FFPDebug(DEBUG_MASS,ffp);

  /*** TKE(A) and its width from the systematics */
  for(int i=0 ; i<ffp->nmass ; i++){
    tke[i] =  FFPSystematics_TKE_A(sf,ffp->zf,ffp->af,i+ffp->mass_first,tke_input);
    dtke[i] = FFPSystematics_TKE_A_Width(sf,ffp->zf,ffp->af,i+ffp->mass_first);
  }

  /*** Z-distribution by Wahl's Zp model */
  WahlZpModel(ffp->zf,ffp->af,ffp->mass_first,ffp->nmass,ffp->ncharge,ffp->ex,charge_dist,charge_first,fzn);

  /*** prepare Z and A list and reaction Q-value in FissionPair object */
  FFPZAList(ffp);
  if(DEBUG_ID & DEBUG_LIST) FFPDebug(DEBUG_LIST,ffp);

  /*** independent yield for a given set of Z and A */
  FFPYield(ffp);
  if(DEBUG_ID & DEBUG_NORM) FFPDebug(DEBUG_NORM,ffp);
  if(DEBUG_ID & DEBUG_CHRG) FFPDebug(DEBUG_CHRG,ffp);

  /*** TXE distributions */
  FFPTXEDistribution(ffp);

  /*** check average TKE and distribution */
  double t0 = tke_input;
  double t1 = FFPAverageTKE(ffp);

  /*** adjust TKE(A) to give the fixed value */
  // if(t0 != 0.0 && t1 != 0.0){
  //   for(int i=0 ; i<ffp->nmass ; i++) tke[i] = FFPSystematics_TKE_A(sf,ffp->zf,ffp->af,i+ffp->mass_first,tke_input) * t0 / t1;
  //   FFPZAList(ffp);
  //   FFPYield(ffp);
  //   FFPTXEDistribution(ffp);
  // }

  if(t0 != 0.0){
    double delta = t0 - t1;
    for(int i=0 ; i<ffp->nmass ; i++) tke[i] = FFPSystematics_TKE_A(sf,ffp->zf,ffp->af,i+ffp->mass_first,tke_input) + delta;
    FFPZAList(ffp);
    FFPYield(ffp);
    FFPTXEDistribution(ffp);
    t1 = FFPAverageTKE(ffp);
  }

  if(DEBUG_ID & DEBUG_ENRG) FFPDebug(DEBUG_ENRG,ffp);
  if(DEBUG_ID & DEBUG_ATKE) FFPDebug(DEBUG_ATKE,ffp);


  /*** free allocated memory */
  FFPDeleteAllocated(ffp);
}


/**********************************************************/
/*      Store Z and A Data                                */
/**********************************************************/
void FFPZAList(FissionFragmentPair *ffp)
{
  double mcn = mass_excess(ffp->zf,ffp->af) + ffp->ex;

  int k = 0;
  for(int i=0 ; i<ffp->nmass ; i++){
    int al = i + ffp->mass_first;
    int ah = ffp->af - al;

    if((double)al > ffp->ac) break;

    for(int j=0 ; j< ffp->ncharge ; j++){
      int zl = j + charge_first[i];
      int zh = ffp->zf - zl;

      double q =  mcn - (mass_excess(zl,al) + mass_excess(zh,ah));

      ffp->fragment[k].setPair(zl,al,zh,ah);
      ffp->fragment[k].qval  = q;
      ffp->fragment[k].yield = 0.0;

      k++;
      if(k > ffp->getMaxPair()){
        cerr << "too many fission pairs " << k << endl;
        return;
      }

      if((al == ah) && (zl == zh)) break;
    }
  }

  ffp->setN(k);
}


/**********************************************************/
/*      Joint Mass and Charge Distributions               */
/**********************************************************/
void FFPYield(FissionFragmentPair *ffp)
{
  /* find min/max Al in FissionPair object */
  int almin = 1000, almax = 0;
  for(int i=0 ; i<ffp->getN() ; i++){
    int al = ffp->fragment[i].getAl();
    if(al < almin) almin = al;
    if(al > almax) almax = al;
  }

  for(int al=almin ; al<=almax ; al++){
    int i = al - ffp->mass_first; // index for Al in the charge_dist array

    /* for each A, determine the elements those are within the Z-range
       and renormalize them */
  
    double sum = 0.0;
    for(int k=0 ; k<ffp->getN() ; k++){
      if(ffp->fragment[k].getAl() != (unsigned int)al) continue;
      int j = ffp->fragment[k].getZl() - charge_first[i]; // index for Zl for a given Al
      if((0 <= j) && (j < ffp->ncharge)) sum += charge_dist[i][j];
    }

    /* mass chain yield times F(Z) */
    if(sum > 0.0) sum = 1.0/sum;
    else sum = 0.0;

    /* when symmetric, divide by 2 to make a smooth mass chain yield */
    if(ffp->af%2 == 0 && al == almax) sum *= 0.5;

    for(int k=0 ; k<ffp->getN() ; k++){
      if(ffp->fragment[k].getAl() != (unsigned int)al) continue;
      int j = ffp->fragment[k].getZl() - charge_first[i];
      if((0 <= j) && (j < ffp->ncharge)) ffp->fragment[k].yield = chain_yield[i] * charge_dist[i][j] * sum;
    }
  }
}


/**********************************************************/
/*      Calculate average TXE                             */
/**********************************************************/
void FFPTXEDistribution(FissionFragmentPair *ffp)
{
  const double TXEcut = 0.5; // eliminate pairs if excitation energy is too low

  for(int k=0 ; k<ffp->getN() ; k++){

    /* average TKE */
    double tke0 = tke[ffp->fragment[k].getAl() - ffp->mass_first];

    /* find the same A and calculate charge-product */
    double s = 0.0;
    int n = 0;
    for(int l=0 ; l<ffp->getN() ; l++){
      if(ffp->fragment[k].getAl() == ffp->fragment[l].getAl()){
        s += (double)ffp->fragment[l].getZl() * (double)ffp->fragment[l].getZh();
        n++;
      }
    }
    /* we assume TKE for each separation is proportional to the Coulomb repulsion */
    s = tke0 * (double)n / s;
    ffp->fragment[k].ek = (double)ffp->fragment[k].getZl() * (double)ffp->fragment[k].getZh() * s;

    /* average TXE */
    ffp->fragment[k].ex = ffp->fragment[k].qval - ffp->fragment[k].ek;

    /* width of TXE is proportional to TKE */
    ffp->fragment[k].sigma =  ffp->fragment[k].ek * dtke[ffp->fragment[k].getAl() - ffp->mass_first] / tke0;

    /* if average TXE is too small, exclude this pair */
    if(ffp->fragment[k].ex < TXEcut) ffp->fragment[k].yield = 0.0;
  }

  /* re-adjust yields */
  double s = 0.0;
  for(int k=0 ; k<ffp->getN() ; k++) s += ffp->fragment[k].yield;
  for(int k=0 ; k<ffp->getN() ; k++) ffp->fragment[k].yield /= s;
}


/**********************************************************/
/*      Calculate Average TKE                             */
/**********************************************************/
double FFPAverageTKE(FissionFragmentPair *ffp)
{
  /* Y sum and average TKE */
  double y = 0.0, e = 0.0;
  for(int k=0 ; k<ffp->getN() ; k++){
    y += ffp->fragment[k].yield;
    e += ffp->fragment[k].yield * ffp->fragment[k].ek;
  }

  return e;
}


/**********************************************************/
/*      Expand List to Include All FFs                    */
/**********************************************************/
void FFPExpandList(const int nc, FissionFragmentPair *ffp)
{
  if(nc == 1) return;

  double mcn = mass_excess(ffp[0].zf,ffp[0].af) + ffp[0].ex;

  /*** ffp[0] will include all the FF pairs */
  
  for(int i=1 ; i<nc ; i++){
    for(int k1=0 ; k1<ffp[i].getN() ; k1++){
      int zl1 = ffp[i].fragment[k1].getZl();
      int al1 = ffp[i].fragment[k1].getAl();
      int zh1 = ffp[i].fragment[k1].getZh();
      int ah1 = ffp[i].fragment[k1].getAh();


      /*** fisrt, search for (Zl,Al) in ffp[0] */
      bool found = false;
      int km = ffp[0].getN();
      for(int k0=0 ; k0<km ; k0++){
        int zl0 = ffp[0].fragment[k0].getZl();
        int al0 = ffp[0].fragment[k0].getAl();
        if( (zl0 == zl1) && (al0 == al1) ){ found = true; break; }
      }

      if( (!found) && (km < ffp[0].getMaxPair()) ){
        int zh0 = ffp[0].zf - zl1;  // this is for the first CN
        int ah0 = ffp[0].af - al1;

        double q =  mcn - (mass_excess(zl1,al1) + mass_excess(zh0,ah0));

        if( (al1 == ah0) && (zl1 > zh0) ) ffp[0].fragment[km].setPair(zh0,al1,zl1,ah0); // special case for Zl > Zh
        else  ffp[0].fragment[km].setPair(zl1,al1,zh0,ah0);
        ffp[0].fragment[km].qval  = q;
        ffp[0].fragment[km].ek    = 0.0;
        ffp[0].fragment[km].ex    = 0.0;
        ffp[0].fragment[km].sigma = 0.0;
        ffp[0].fragment[km].yield = 0.0;
        km ++;
        ffp[0].setN(km);
      }

      /*** second, search for (Zh,Ah) in ffp[0] */
      found = false;
      for(int k0=0 ; k0<km ; k0++){
        int zh0 = ffp[0].fragment[k0].getZh();
        int ah0 = ffp[0].fragment[k0].getAh();
        if( (zh0 == zh1) && (ah0 == ah1) ){ found = true; break; }
      }

      if( (!found) && (km < ffp[0].getMaxPair()) ){
        int zl0 = ffp[0].zf - zh1;
        int al0 = ffp[0].af - ah1;

        double q =  mcn - (mass_excess(zl0,al0) + mass_excess(zh1,ah1));

        if( (al0 == ah1) && (zl0 > zh1) ) ffp[0].fragment[km].setPair(zh1,al0,zl0,ah1);
        else ffp[0].fragment[km].setPair(zl0,al0,zh1,ah1);
        ffp[0].fragment[km].qval  = q;
        ffp[0].fragment[km].ek    = 0.0;
        ffp[0].fragment[km].ex    = 0.0;
        ffp[0].fragment[km].sigma = 0.0;
        ffp[0].fragment[km].yield = 0.0;
        km ++;
        ffp[0].setN(km);
      }
    }
  }
//FFPDebug(DEBUG_NORM,ffp);

  /*** duplication test */
/*
  for(int k0=0 ; k0<ffp[0].getN() ; k0++){
    int zl0 = ffp[0].fragment[k0].getZl();
    int al0 = ffp[0].fragment[k0].getAl();
    int zh0 = ffp[0].fragment[k0].getZh();
    int ah0 = ffp[0].fragment[k0].getAh();

    for(int k1=k0+1 ; k1<ffp[0].getN() ; k1++){
      int zl1 = ffp[0].fragment[k1].getZl();
      int al1 = ffp[0].fragment[k1].getAl();
      int zh1 = ffp[0].fragment[k1].getZh();
      int ah1 = ffp[0].fragment[k1].getAh();

      if( (zl0 == zl1) && (al0 == al1) && (zh0 == zh1) && (ah0 == ah1) ){
        cout << "duplicated " << zl0 << " " << al0 << " " << zh0 << " " << ah0 << endl;
      }
    }
  }
*/
}


/**********************************************************/
/*      Allocate Temporal Arrays                          */
/**********************************************************/
void FFPAllocateMemory(FissionFragmentPair *ffp)
{
  try{
    tke          = new double [ffp->nmass];   // TKE(A)
    dtke         = new double [ffp->nmass];   // width of TKE(A)
    chain_yield  = new double [ffp->nmass];   // Chain Yield Y(A)
    charge_first = new int    [ffp->nmass];   // lowest Z number for each Z-dist
    charge_dist  = new double * [ffp->nmass]; // probability of Z for a given A
    for(int i=0 ; i<ffp->nmass; i++) charge_dist[i] = new double [ffp->ncharge];
  }
  catch(bad_alloc &e){
    message << "memory allocation error";
    cohTerminateCode("FFPAllocateMemory");
  }

  for(int a=0 ; a<ffp->nmass; a++){
    tke[a] = dtke[a] = chain_yield[a] = 0.0;
    charge_first[a] = 0;
    for(int z=0 ; z<ffp->ncharge; z++) charge_dist[a][z] = 0.0;
  }
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void FFPDeleteAllocated(FissionFragmentPair *ffp)
{
  delete [] tke;
  delete [] dtke;
  delete [] chain_yield;

  for(int a=0 ; a<ffp->nmass; a++) delete [] charge_dist[a];
  delete [] charge_dist;
  delete [] charge_first;
}


/**********************************************************/
/*      for Debugging                                     */
/**********************************************************/
void FFPDebug(unsigned char id, FissionFragmentPair *ffp)
{
  int zlmin = 1000, zlmax = 0;
  int almin = 1000, almax = 0;
  for(int i=0 ; i<ffp->getN() ; i++){
    int zl = ffp->fragment[i].getZl();
    if(zl < zlmin) zlmin = zl;
    if(zl > zlmax) zlmax = zl;
    int al = ffp->fragment[i].getAl();
    if(al < almin) almin = al;
    if(al > almax) almax = al;
  }


  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  /* print mass chain yield */
  if(id & DEBUG_MASS){
    double s = 0.0;
    for(int i=0 ; i<ffp->nmass ; i++){
      s += chain_yield[i];
      cout << setw(5) << i + ffp->mass_first;
      cout << setw(14) << chain_yield[i];
      cout << setw(14) << s << endl;
    }
  }

  /* print all list */
  if(id & DEBUG_LIST){
    for(int k=0 ; k<ffp->getN() ; k++){
      cout << setw(5) << k;
      cout << setw(5) << ffp->fragment[k].getZl();
      cout << setw(5) << ffp->fragment[k].getAl();
      cout << setw(5) << ffp->fragment[k].getZh();
      cout << setw(5) << ffp->fragment[k].getAh();
      cout << setw(14) << ffp->fragment[k].qval << endl;
    }
  }

  /* check normalization */
  if(id & DEBUG_NORM){
    double s = 0.0;
    for(int k=0 ; k<ffp->getN() ; k++){
      s += ffp->fragment[k].yield;
      cout << setw(5) << k;
      cout << setw(5) << ffp->fragment[k].getZl();
      cout << setw(5) << ffp->fragment[k].getAl();
      cout << setw(5) << ffp->fragment[k].getZh();
      cout << setw(5) << ffp->fragment[k].getAh();
      cout << setw(14) << ffp->fragment[k].yield;
      cout << setw(14) << s << endl;
    }
  }

  /* check charge distribution */
  if(id & DEBUG_CHRG){
    double x = 0.0;
    for(int zl=zlmin ; zl<=zlmax ; zl++){
      double y = 0.0;
      for(int k=0 ; k<ffp->getN() ; k++){
        if(ffp->fragment[k].getZl() == (unsigned int)zl) y += ffp->fragment[k].yield;
      }
      x += y;
      cout << setw(5) << zl;
      cout << setw(5) << ffp->zf - zl;
      cout << setw(14) << y;
      cout << setw(14) << x << endl;
    }
  }

  /* check TKE and TXE */
  if(id & DEBUG_ENRG){
    double s = 0.0;
    for(int k=0 ; k<ffp->getN() ; k++){
      if(ffp->fragment[k].yield == 0.0) continue;
      s += ffp->fragment[k].yield;
      cout << setw(5)  << k;
      cout << setw(5)  << ffp->fragment[k].getZl();
      cout << setw(5)  << ffp->fragment[k].getAl();
      cout << setw(5)  << ffp->fragment[k].getZh();
      cout << setw(5)  << ffp->fragment[k].getAh();
      cout << setw(14) << ffp->fragment[k].qval;
      cout << setw(14) << ffp->fragment[k].ek;
      cout << setw(14) << ffp->fragment[k].ex;
      cout << setw(14) << ffp->fragment[k].yield;
      cout << setw(14) << s << endl;
    }
  }

  /* check TKE and its width */
  if(id & DEBUG_ATKE){
    for(int al=almin ; al<=almax ; al++){
      double tke0 = tke[al - ffp->mass_first];

      double t = 0.0, w = 0.0, y = 0.0;
      for(int k=0 ; k<ffp->getN() ; k++){
        if(ffp->fragment[k].getAl() == (unsigned int)al){
          t += ffp->fragment[k].ek * ffp->fragment[k].yield;
          w += ffp->fragment[k].sigma * ffp->fragment[k].sigma * ffp->fragment[k].yield;
          y += ffp->fragment[k].yield;
        }
      }
      if(y > 0.0){
        t = t/y;
        w = sqrt(w/y);
      }
      cout << setw(5)  << al;
      cout << setw(5)  << ffp->af - al;
      cout << setw(14) << tke0;
      cout << setw(14) << t;
      cout << setw(14) << t/tke0;
      cout << setw(14) << w << endl;
    }
  }
}
