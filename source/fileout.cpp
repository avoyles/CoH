/******************************************************************************/
/*  fileout.cpp                                                               */
/*        write calculated results on files                                   */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

#include "structur.h"
#include "nucleus.h"
#include "fileout.h"
#include "global.h"
#include "terminate.h"

#undef LONGOUTPUT

static const int N_OUTPUT_CROSSSEC =  6;
static const int N_OUTPUT_SIGPART  = 24;
static const int N_OUTPUT_FISSION  =  5;
static const int N_PART_KIND       =  7;
static const int N_OUTPUT_LEVELS   = 42;

#ifdef LONGOUTPUT
static const int COLUMN_WIDTH      = 21;
static const int DECIMAL_WIDTH     = 14;
#else
static const int COLUMN_WIDTH      = 13;
static const int DECIMAL_WIDTH     =  6;
#endif

static void    fioWriteCrossSection1           (const bool, const string, int, double);
static void    fioWriteCrossSection2           (const bool, const string, int, double);
static void    fioWriteCrossSection3           (const bool, const string, int, double);
static void    fioWriteLevelExcite             (const bool, const string, double *, double, int, Nucleus *);
static void    fioWriteAngularDistribution     (const string, double *, double);
static void    fioWriteLegendreCoefficient     (const string, double *, double, int);

static inline void  fioSetLevelEnergy          (double *, Nucleus *);
static inline bool  fioCheckExist              (string);

static string rname1[N_OUTPUT_CROSSSEC] = {
  " Total       "," Elastic     "," Non-elastic "," Capture     "," Inelastic   "," Fission     "};

static string rname2[N_OUTPUT_SIGPART ] = {
  " (x,n)       "," (x,p)       "," (x,A)       "," (x,d)       "," (x,t)       "," (x,h)       ",
  " (x,2n)      "," (x,np)      "," (x,nA)      "," (x,nd)      "," (x,nt)      "," (x,nh)      ",
  " (x,2p)      "," (x,pA)      "," (x,pd)      "," (x,2A)      "," (x,Ad)      "
  " (x,3n)      "," (x,2np)     "," (x,2nA)     "," (x,npA)     "," (x,n2p)     ",
  " (x,4n)      "," (x,5n)      "};

static string rname3[N_OUTPUT_FISSION ] = {
  " (x,0nf)     "," (x,1nf)     "," (x,2nf)     "," (x,3nf)     "," (x,4nf)     "};

  

/*** lower limit of output values, enforce zero */
static const double eps = 1.0e-99;

/**********************************************************/
/*     Output Cross Sections and Angular Distributions    */
/**********************************************************/
void cohFileOut(int nc, int tid, double elab)
// nc    : number of compound nucleus
// tid   : index for the target
// elab  : incident LAB energy 
{
  string fname;
  ostringstream fn;
  bool fexist = false;
  double elev[N_OUTPUT_LEVELS];

  fioSetLevelEnergy(elev,&ncl[tid]);

  if(prn.angdist){
    fn.str("");
    fn << filename_angular_distribution << "." << file_extension;
    fioWriteAngularDistribution(fn.str(), elev, elab);
  }
  else{
    /*** check if cross section output files exist, if not,  */
    fn.str("");
    fn << filename_cross_section << "." << file_extension;
    fexist = fioCheckExist(fn.str());
    fioWriteCrossSection1(fexist, fn.str(), nc, elab);

    fn.str("");
    fn << filename_particle_production << "." << file_extension;
    fexist = fioCheckExist(fn.str());
    fioWriteCrossSection2(fexist, fn.str(), nc, elab);

    if(ctl.fission){
      fn.str("");
      fn << filename_fission << "." << file_extension;
      fexist = fioCheckExist(fn.str());
      fioWriteCrossSection3(fexist, fn.str(), nc, elab);
    }
  }

  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;

    int p = ncl[0].cdt[id].next;
    if(ncl[p].ndisc == 0) continue;

    fioSetLevelEnergy(elev,&ncl[p]);

    fn.str("");
    if(prn.angdist){
      fn << filename_legendre_coefficient << setw(1) << id << "." << file_extension;
      fioWriteLegendreCoefficient(fn.str(), elev, elab, id);
    }
    else{
      fn << filename_level_excite << setw(1) << id << "." << file_extension;
      fexist = fioCheckExist(fn.str());
      fioWriteLevelExcite(fexist, fn.str(), elev, elab, id, &ncl[p]);
    }
  }
}


/**********************************************************/
/*     Store Discrete Level Energies in Array             */
/**********************************************************/
void fioSetLevelEnergy(double *elev, Nucleus *n)
{
  for(int k=0 ; k<N_OUTPUT_LEVELS ; k++) elev[k] = 0.0;

  int m = N_OUTPUT_LEVELS-1;
  if(n->ndisc < m) m = n->ndisc;

  for(int k=1 ; k<m ; k++) elev[k] = n->lev[k].energy;
}


/**********************************************************/
/*     Check If File Already Exists in the Current Dir    */
/**********************************************************/
bool fioCheckExist(string fname)
{
  ifstream fin;
  bool     fexist = false;

  fin.open(&fname[0],ios::in);
  if(!fin) fexist = false;
  else     fexist = true;
  fin.close();

  return(fexist);
}


/**********************************************************/
/*     Output Selected Cross Sections on File             */
/**********************************************************/
void fioWriteCrossSection1(const bool fexist, const string fname, int nc, double e)
{
  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSection1");
  }

  char   s[N_PART_KIND+1] = "0000000";
  string str;
  double c1[N_OUTPUT_CROSSSEC];

  for(int i=0 ; i<N_OUTPUT_CROSSSEC ; i++) c1[i]=0.0;

  c1[0] = crx.total;                    // total cross section
  c1[1] = crx.elastic + crx.compound;   // elastic scattering cross section
  c1[2] = crx.reaction;                 // non-elastic cross section

  /*** capture, inelastic, and fission */
  for(int i=0 ; i<nc ; i++) {
    for(int j=0 ; j<min(N_PART_KIND,MAX_CHANNEL) ; j++)
      sprintf(&s[j],"%1d",(int)crx.prod[i].par[j]);
    str = (string)s;

    if(     str=="0000000"){
      c1[3] = crx.prod[i].xsec;
    }
    else if(str=="0100000"){
      double z =crx.prod[i].xsec - crx.compound;
      c1[4] = (z < 1.0e-10) ? 0.0 : z;
    }
    c1[5] += crx.prod[i].fiss;
  }

  if(!fexist){
    fp << "#" << setw(COLUMN_WIDTH-1) << "            ";
    for(int i=0 ; i<N_OUTPUT_CROSSSEC ; i++) fp << setw(COLUMN_WIDTH) << rname1[i];
    fp << endl;
  }

  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  fp << setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_CROSSSEC ; i++){
    fp << setw(COLUMN_WIDTH) << ( (c1[i] < eps) ? 0.0 : c1[i] );
  }
  fp << endl;
  fp.close();
}


/**********************************************************/
/*     Output Particle Production Cross Sections          */
/**********************************************************/
void fioWriteCrossSection2(const bool fexist, const string fname, int nc, double e)
{
  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSection2");
  }

  char   s[N_PART_KIND+1] = "0000000";
  string str;
  double c2[N_OUTPUT_SIGPART];

  for(int i=0 ; i<N_OUTPUT_SIGPART ; i++) c2[i]=0.0;

  for(int i=0 ; i<nc ; i++) {
    for(int j=0 ; j<min(N_PART_KIND,MAX_CHANNEL) ; j++)
      sprintf(&s[j],"%1x",(int)crx.prod[i].par[j]);
    str = (string)s;

    if(     str=="0100000"){ c2[ 0] = crx.prod[i].xsec; }  // (x,n)
    else if(str=="0010000"){ c2[ 1] = crx.prod[i].xsec; }  // (x,p)
    else if(str=="0001000"){ c2[ 2] = crx.prod[i].xsec; }  // (x,A)
    else if(str=="0000100"){ c2[ 3] = crx.prod[i].xsec; }  // (x,d)
    else if(str=="0000010"){ c2[ 4] = crx.prod[i].xsec; }  // (x,t)
    else if(str=="0000001"){ c2[ 5] = crx.prod[i].xsec; }  // (x,h)

    else if(str=="0200000"){ c2[ 6] = crx.prod[i].xsec; }  // (x,2n)
    else if(str=="0110000"){ c2[ 7] = crx.prod[i].xsec; }  // (x,np)
    else if(str=="0101000"){ c2[ 8] = crx.prod[i].xsec; }  // (x,nA)
    else if(str=="0100100"){ c2[ 9] = crx.prod[i].xsec; }  // (x,nd)
    else if(str=="0100010"){ c2[10] = crx.prod[i].xsec; }  // (x,nt)
    else if(str=="0100001"){ c2[11] = crx.prod[i].xsec; }  // (x,nh)

    else if(str=="0020000"){ c2[12] = crx.prod[i].xsec; }  // (x,2p)
    else if(str=="0011000"){ c2[13] = crx.prod[i].xsec; }  // (x,pA)
    else if(str=="0010100"){ c2[14] = crx.prod[i].xsec; }  // (x,pd)
    else if(str=="0002000"){ c2[15] = crx.prod[i].xsec; }  // (x,2A)
    else if(str=="0001100"){ c2[16] = crx.prod[i].xsec; }  // (x,Ad)

    else if(str=="0300000"){ c2[17] = crx.prod[i].xsec; }  // (x,3n)
    else if(str=="0210000"){ c2[18] = crx.prod[i].xsec; }  // (x,2np)
    else if(str=="0201000"){ c2[19] = crx.prod[i].xsec; }  // (x,2nA)
    else if(str=="0111000"){ c2[20] = crx.prod[i].xsec; }  // (x,npA)
    else if(str=="0120000"){ c2[21] = crx.prod[i].xsec; }  // (x,n2p)

    else if(str=="0400000"){ c2[22] = crx.prod[i].xsec; }  // (x,4n)
    else if(str=="0500000"){ c2[23] = crx.prod[i].xsec; }  // (x,5n)
  }

  if(!fexist){
    fp << "#" << setw(COLUMN_WIDTH-1) << "            ";
    for(int i=0 ; i<N_OUTPUT_SIGPART ; i++) fp << setw(COLUMN_WIDTH) << rname2[i];
    fp << endl;
  }

  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  fp << setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_SIGPART ; i++){
    fp << setw(COLUMN_WIDTH) << ( (c2[i] < eps) ? 0.0 : c2[i] );
  }
  fp << endl;
  fp.close();
}


/**********************************************************/
/*     Output Fission Cross Sections on File              */
/**********************************************************/
void fioWriteCrossSection3(const bool fexist, const string fname, int nc, double e)
{
  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSection");
  }

  double c3[N_OUTPUT_FISSION];
  
  for(int i=0 ; i<N_OUTPUT_FISSION ; i++) c3[i]=0.0;
  for(int i=0 ; i<nc ; i++){
    if(i >= N_OUTPUT_FISSION) break;
    c3[i] = crx.prod[i].fiss;
  }

  if(!fexist){
    fp << "#" << setw(COLUMN_WIDTH-1) << "            ";
    for(int i=0 ; i<N_OUTPUT_FISSION ; i++) fp << setw(COLUMN_WIDTH) << rname3[i];
    fp << endl;
  }

  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  fp << setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_FISSION ; i++){
    fp << setw(COLUMN_WIDTH) << ( (c3[i] < eps) ? 0.0 : c3[i] );
  }
  fp << endl;
  fp.close();
}


/**********************************************************/
/*     Output Level Excitation Cross Sections on File     */
/**********************************************************/
void fioWriteLevelExcite(const bool fexist, const string fname, double *y, double e, int id, Nucleus *n)
{
  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteInelastic");
  }

  double x[N_OUTPUT_LEVELS];

  for(int k=0 ; k<N_OUTPUT_LEVELS ; k++) x[k] = 0.0;
  
  /*** sum all level inelastic up to min(N_OUTPUT_LEVELS-1,ndisc) */
  int m = N_OUTPUT_LEVELS-1;
  if(n->ndisc < m) m = n->ndisc;

  double sum = 0.0;
  for(int k=0 ; k<m ; k++){
    sum += (x[k] = crx.levexcite[id][k]);
  }

  /*** the last element contains continuum inelastic cross section */
  x[N_OUTPUT_LEVELS-1] = n->lpop[0] - sum;
  y[N_OUTPUT_LEVELS-1] = y[m-1] + 0.001;

  /*** not so exact for the continuum, to avoid round off error    */
  if(x[N_OUTPUT_LEVELS-1] < 0.0) x[N_OUTPUT_LEVELS-1] = 0.0;

  /*** if the c/s to the highest level is zero, no continuum can be assumed */
  if(x[m-1] == 0.0) x[N_OUTPUT_LEVELS-1] = 0.0;

  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  if(!fexist){
    fp << "#" << setw(COLUMN_WIDTH-1) << "            ";
    for(int i=0 ; i<N_OUTPUT_LEVELS ; i++) fp << setw(COLUMN_WIDTH) << y[i];
    fp << endl;
  }
  fp << setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_LEVELS ; i++) fp << setw(COLUMN_WIDTH) << x[i];
  fp << endl;
  fp.close();
}


/**********************************************************/
/*     Output Scattering Angular Distributions on File    */
/**********************************************************/
void fioWriteAngularDistribution(string fname, double *y, double e)
{
  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteAngularDistribution");
  }
  
  fp << setprecision(DECIMAL_WIDTH-1);
  fp << setiosflags(ios::scientific);
  fp << "#" << setw(COLUMN_WIDTH-1) << e << setw(COLUMN_WIDTH) << "   Angdist   ";

  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  for(int k=1 ; k<N_OUTPUT_LEVELS-1 ; k++) fp << setw(COLUMN_WIDTH) << y[k];
  fp << endl;

  for(int i=0 ; i<MAX_ANGDIST ; i++){
    fp << setw(COLUMN_WIDTH) << crx.theta[i];
    for(int j=0 ; j<N_OUTPUT_LEVELS-1 ; j++) fp << setw(COLUMN_WIDTH) << crx.angdist[j][i];
    fp << endl;
  }

  fp.close();
}


/**********************************************************/
/*     Output Legendre Coefficients on File               */
/**********************************************************/
void fioWriteLegendreCoefficient(string fname, double *y, double e, int id)
{
  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteLegendreCoefficient");
  }

  fp << setprecision(DECIMAL_WIDTH-1);
  fp << setiosflags(ios::scientific);
  fp << "#" << setw(COLUMN_WIDTH-1) << e << setw(COLUMN_WIDTH) << "   Legcoef   ";

  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  for(int k=1 ; k<N_OUTPUT_LEVELS-1 ; k++) fp << setw(COLUMN_WIDTH) << y[k];
  fp << endl;

  fp << setprecision(DECIMAL_WIDTH-1);
  for(int i=0 ; i<MAX_J ; i++){
    fp << setw(COLUMN_WIDTH) << i;
    for(int j=0 ; j<N_OUTPUT_LEVELS-1 ; j++) fp << setw(COLUMN_WIDTH) <<  crx.legcoef[id][j][i];
    fp << endl;
  }

  fp.close();
}
