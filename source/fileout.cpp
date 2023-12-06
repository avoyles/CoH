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
#include "metastable.h"
#include "terminate.h"

static const int N_PART_KIND       =  7;
static const int N_OUTPUT_LEVELS   = 42;

#undef LONGOUTPUT
#ifdef LONGOUTPUT
static const int COLUMN_WIDTH      = 21;
static const int DECIMAL_WIDTH     = 14;
#else
static const int COLUMN_WIDTH      = 13;
static const int DECIMAL_WIDTH     =  6;
#endif

static void    fioWriteCrossSectionMain        (const bool, const string, const int, const double);
static void    fioWriteCrossSectionParticle    (const bool, const string, const int, const double);
static void    fioWriteCrossSectionIsomer      (const bool, const string, const int, const double);
static void    fioWriteCrossSectionFission     (const bool, const string, const int, const double);
static void    fioWriteLevelExcite             (const bool, const string, double *, const double, const int, Nucleus *);
static void    fioWriteAngularDistribution     (const string, double *, const double);
static void    fioWriteLegendreCoefficient     (const string, double *, const double, const int);

static inline void  fioSetLevelEnergy          (double *, Nucleus *);
static inline bool  fioCheckExist              (string);


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
  double elev[N_OUTPUT_LEVELS];

  /*** print out angular distributions for elastic and inelastic scattering */
  if(prn.angdist){

    fioSetLevelEnergy(elev,&ncl[tid]);

    fn.str("");
    fn << filename_angular_distribution << "." << file_extension;
    fioWriteAngularDistribution(fn.str(), elev, elab);
  }

  /*** print out cross sections */
  else{
    /*** main cross sections, total, elastic, capture, inelasitc, fission  */
    fn.str("");
    fn << filename_cross_section << "." << file_extension;
    fioWriteCrossSectionMain(fioCheckExist(fn.str()), fn.str(), nc, elab);

    /*** all reaction cross sections */
    fn.str("");
    fn << filename_particle_production << "." << file_extension;
    fioWriteCrossSectionParticle(fioCheckExist(fn.str()), fn.str(), nc, elab);

    /*** isomeric ratios */
    fn.str("");
    fn << filename_isomeric_ratio << "." << file_extension;
    fioWriteCrossSectionIsomer(fioCheckExist(fn.str()), fn.str(), nc, elab);

    /*** multi-chance fission cross sections */
    if(ctl.fission){
      fn.str("");
      fn << filename_fission << "." << file_extension;
      fioWriteCrossSectionFission(fioCheckExist(fn.str()), fn.str(), nc, elab);
    }
  }

  /*** for binary reactions, scattering to discrete levels */
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
      fioWriteLevelExcite(fioCheckExist(fn.str()), fn.str(), elev, elab, id, &ncl[p]);
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
void fioWriteCrossSectionMain(const bool fexist, const string fname, const int nc, const double e)
{
  const int N_OUTPUT_CROSSSEC = 6;
  static const string rname[N_OUTPUT_CROSSSEC] = {
  " Total       "," Elastic     "," CNFormation "," Capture     "," Inelastic   "," Fission     "};

  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSectionMain");
  }

  char   s[N_PART_KIND+1] = "0000000";
  string str;
  double cx[N_OUTPUT_CROSSSEC];

  for(int i=0 ; i<N_OUTPUT_CROSSSEC ; i++) cx[i]=0.0;

  cx[0] = crx.total;                    // total cross section
  cx[1] = crx.elastic + crx.compound;   // elastic scattering cross section
  cx[2] = crx.reaction;                 // compond formation cross section

  /*** capture, inelastic, and fission */
  for(int i=0 ; i<nc ; i++) {
    for(int j=0 ; j<min(N_PART_KIND,MAX_CHANNEL) ; j++)
      sprintf(&s[j],"%1d",(int)crx.prod[i].par[j]);
    str = (string)s;

    if(     str=="0000000"){
      cx[3] = crx.prod[i].xsec;
    }
    else if(str=="0100000"){
      double z =crx.prod[i].xsec - crx.compound;
      cx[4] = (z < 1.0e-10) ? 0.0 : z;
    }
    cx[5] += crx.prod[i].fiss;
  }

  /*** print header if file not exist */
  if(!fexist){
    fp << "#" << setw(COLUMN_WIDTH-1) << " Reaction    ";
    for(int i=0 ; i<N_OUTPUT_CROSSSEC ; i++) fp << setw(COLUMN_WIDTH) << rname[i];
    fp << endl;
  }

  /*** print cross sections */
  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  fp << setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_CROSSSEC ; i++){
    fp << setw(COLUMN_WIDTH) << ( (cx[i] < eps) ? 0.0 : cx[i] );
  }
  fp << endl;
  fp.close();
}


/**********************************************************/
/*     Output Particle Production Cross Sections          */
/**********************************************************/
void fioWriteCrossSectionParticle(const bool fexist, const string fname, const int nc, const double e)
{
  const int N_OUTPUT_SIGPART = 24;
  static const string rname[N_OUTPUT_SIGPART ] = {
  " (x,n)       "," (x,p)       "," (x,A)       "," (x,d)       "," (x,t)       "," (x,h)       ",
  " (x,2n)      "," (x,np)      "," (x,nA)      "," (x,nd)      "," (x,nt)      "," (x,nh)      ",
  " (x,2p)      "," (x,pA)      "," (x,pd)      "," (x,2A)      "," (x,Ad)      "
  " (x,3n)      "," (x,2np)     "," (x,2nA)     "," (x,npA)     "," (x,n2p)     ",
  " (x,4n)      "," (x,5n)      "};

  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSectionParticle");
  }

  char   s[N_PART_KIND+1] = "0000000";
  string str;
  double cx[N_OUTPUT_SIGPART];

  for(int i=0 ; i<N_OUTPUT_SIGPART ; i++) cx[i] = 0.0;

  for(int i=0 ; i<nc ; i++) {
    for(int j=0 ; j<min(N_PART_KIND,MAX_CHANNEL) ; j++) sprintf(&s[j],"%1x",(int)crx.prod[i].par[j]);
    str = (string)s;

    if(     str == "0100000"){ cx[ 0] = crx.prod[i].xsec; }  // (x,n)
    else if(str == "0010000"){ cx[ 1] = crx.prod[i].xsec; }  // (x,p)
    else if(str == "0001000"){ cx[ 2] = crx.prod[i].xsec; }  // (x,A)
    else if(str == "0000100"){ cx[ 3] = crx.prod[i].xsec; }  // (x,d)
    else if(str == "0000010"){ cx[ 4] = crx.prod[i].xsec; }  // (x,t)
    else if(str == "0000001"){ cx[ 5] = crx.prod[i].xsec; }  // (x,h)

    else if(str == "0200000"){ cx[ 6] = crx.prod[i].xsec; }  // (x,2n)
    else if(str == "0110000"){ cx[ 7] = crx.prod[i].xsec; }  // (x,np)
    else if(str == "0101000"){ cx[ 8] = crx.prod[i].xsec; }  // (x,nA)
    else if(str == "0100100"){ cx[ 9] = crx.prod[i].xsec; }  // (x,nd)
    else if(str == "0100010"){ cx[10] = crx.prod[i].xsec; }  // (x,nt)
    else if(str == "0100001"){ cx[11] = crx.prod[i].xsec; }  // (x,nh)

    else if(str == "0020000"){ cx[12] = crx.prod[i].xsec; }  // (x,2p)
    else if(str == "0011000"){ cx[13] = crx.prod[i].xsec; }  // (x,pA)
    else if(str == "0010100"){ cx[14] = crx.prod[i].xsec; }  // (x,pd)
    else if(str == "0002000"){ cx[15] = crx.prod[i].xsec; }  // (x,2A)
    else if(str == "0001100"){ cx[16] = crx.prod[i].xsec; }  // (x,Ad)

    else if(str == "0300000"){ cx[17] = crx.prod[i].xsec; }  // (x,3n)
    else if(str == "0210000"){ cx[18] = crx.prod[i].xsec; }  // (x,2np)
    else if(str == "0201000"){ cx[19] = crx.prod[i].xsec; }  // (x,2nA)
    else if(str == "0111000"){ cx[20] = crx.prod[i].xsec; }  // (x,npA)
    else if(str == "0120000"){ cx[21] = crx.prod[i].xsec; }  // (x,n2p)

    else if(str == "0400000"){ cx[22] = crx.prod[i].xsec; }  // (x,4n)
    else if(str == "0500000"){ cx[23] = crx.prod[i].xsec; }  // (x,5n)
  }

  /*** print header if file not exist */
  if(!fexist){
    fp << "#" << setw(COLUMN_WIDTH-1) << " Reaction   ";
    for(int i=0 ; i<N_OUTPUT_SIGPART ; i++) fp << setw(COLUMN_WIDTH) << rname[i];
    fp << endl;
  }

  /*** print cross sections */
  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  fp << setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_SIGPART ; i++){
    fp << setw(COLUMN_WIDTH) << ( (cx[i] < eps) ? 0.0 : cx[i] );
  }
  fp << endl;
  fp.close();
}


/**********************************************************/
/*     Output Isomeric State Production Cross Sections    */
/**********************************************************/
void fioWriteCrossSectionIsomer(const bool fexist, const string fname, const int nc, const double e)
{
  const int metamax = 2; // number of isomeric states to be printed
  const int N_OUTPUT_ISOMER = 13;
  static const string name[N_OUTPUT_ISOMER] = {
      "0000000", "0100000", "0010000", "0001000", "0000100", "0000010", "0000001",
      "0200000", "0110000", "0101000", "0100100", "0100010", "0100001"};

  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSectionIsomer");
  }

  char   s[N_PART_KIND+1] = "0000000";
  string str;
  int    meta[N_OUTPUT_ISOMER];

  double **cm = new double * [N_OUTPUT_ISOMER];
  double **em = new double * [N_OUTPUT_ISOMER];
  double **tm = new double * [N_OUTPUT_ISOMER];
  for(int m=0 ; m<N_OUTPUT_ISOMER ; m++){
    cm[m] = new double [metamax];
    em[m] = new double [metamax];
    tm[m] = new double [metamax];
    meta[m] = 0;
  }

  /*** for each residual, see the reaction type */
  for(int i=0 ; i<nc ; i++) {
    for(int j=0 ; j<min(N_PART_KIND,MAX_CHANNEL) ; j++) sprintf(&s[j],"%1x",(int)crx.prod[i].par[j]);
    str = (string)s;

    /*** check if this reaction will be printed (listed in name[]) */
    int m = -1;
    for(int j=0 ; j<N_OUTPUT_ISOMER ; j++){
       if(str == name[j]){ m = j; break; }
    }
    if(m < 0) continue;

    /*** extract isomeric state production ratios */
    for(int k=1 ; k<=ncl[i].ndisc ; k++){
      if( (ncl[i].lev[k].halflife > thalfmin) && (meta[m] < metamax) ){
        /*** excitation energy of the state */
        em[m][meta[m]] = ncl[i].lev[k].energy;

        /*** half-life of the state */
        tm[m][meta[m]] = ncl[i].lev[k].halflife;

        /*** production cross section ratio to g.s. */
        cm[m][meta[m]] = (ncl[i].lpop[0] > 0.0) ? ncl[i].lpop[k] / ncl[i].lpop[0] : 0.0;

        /* increment meta-stable state counter */
        meta[m]++;
      }
    }
  }

  /*** print header if file not exist */
  if(!fexist){
    fp << "#" << setw(COLUMN_WIDTH-1) << " PartEmit   ";
    for(int i=0 ; i<N_OUTPUT_ISOMER ; i++){
      for(int m=0 ; m<metamax ; m++){
        fp << setw(COLUMN_WIDTH-2) << name[i] << '-' << setw(1) << m+1;
      }
    }
    fp << endl;

    fp << setprecision(DECIMAL_WIDTH);
    fp << setiosflags(ios::scientific);
    fp << "#" << setw(COLUMN_WIDTH-1) << " Excitation ";
    for(int i=0 ; i<N_OUTPUT_ISOMER ; i++){
      for(int m=0 ; m<metamax ; m++){
        fp << setw(COLUMN_WIDTH) << em[i][m];
      }
    }
    fp << endl;

    fp << setprecision(DECIMAL_WIDTH);
    fp << setiosflags(ios::scientific);
    fp << "#" << setw(COLUMN_WIDTH-1) << " HalfLife   ";
    for(int i=0 ; i<N_OUTPUT_ISOMER ; i++){
      for(int m=0 ; m<metamax ; m++){
        fp << setw(COLUMN_WIDTH) << tm[i][m];
      }
    }
    fp << endl;
  }

  /*** print cross sections */
  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  fp << setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_ISOMER ; i++){
    for(int m=0 ; m<metamax ; m++){
      fp << setw(COLUMN_WIDTH) << ( (cm[i][m] < eps) ? 0.0 : cm[i][m] );
    }
  }
  fp << endl;
  fp.close();

  for(int i=0 ; i<N_OUTPUT_ISOMER ; i++){
    delete [] cm[i];
    delete [] em[i];
    delete [] tm[i];
  }
  delete [] cm;
  delete [] em;
  delete [] tm;
}


/**********************************************************/
/*     Output Fission Cross Sections on File              */
/**********************************************************/
void fioWriteCrossSectionFission(const bool fexist, const string fname, const int nc, const double e)
{
  const int N_OUTPUT_FISSION = 5;
  static const string rname[N_OUTPUT_FISSION] = {
  " (x,0nf)     "," (x,1nf)     "," (x,2nf)     "," (x,3nf)     "," (x,4nf)     "};

  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSectionFission");
  }

  double cx[N_OUTPUT_FISSION];
  
  for(int i=0 ; i<N_OUTPUT_FISSION ; i++) cx[i] = 0.0;
  for(int i=0 ; i<nc ; i++){
    if(i >= N_OUTPUT_FISSION) break;
    cx[i] = crx.prod[i].fiss;
  }

  /*** print header if file not exist */
  if(!fexist){
    fp << "#" << setw(COLUMN_WIDTH-1) << " Fission    ";
    for(int i=0 ; i<N_OUTPUT_FISSION ; i++) fp << setw(COLUMN_WIDTH) << rname[i];
    fp << endl;
  }

  /*** print cross sections */
  fp << setprecision(DECIMAL_WIDTH);
  fp << setiosflags(ios::scientific);
  fp << setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_FISSION ; i++){
    fp << setw(COLUMN_WIDTH) << ( (cx[i] < eps) ? 0.0 : cx[i] );
  }
  fp << endl;
  fp.close();
}


/**********************************************************/
/*     Output Level Excitation Cross Sections on File     */
/**********************************************************/
void fioWriteLevelExcite(const bool fexist, const string fname, double *y, const double e, const int id, Nucleus *n)
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
void fioWriteAngularDistribution(string fname, double *y, const double e)
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
void fioWriteLegendreCoefficient(string fname, double *y, const double e, const int id)
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
