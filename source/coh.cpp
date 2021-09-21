/******************************************************************************/
/**                                                                          **/
/**     C o H 3 : The Hauser-Feshbach Code                                   **/
/**                                              Version 3.5 Miranda  (2015) **/
/**                                              T.Kawano                    **/
/**     History                                                              **/
/**       3.0 Callisto : developing version for full Hauser-Feshbach  (2009) **/
/**       3.1 Ariel    : fission modeling                             (2010) **/
/**       3.2 Umbriel  : COH + ECLIPSE unified version                (2012) **/
/**       3.3 Titania  : advanced memory management version           (2013) **/
/**       3.4 Oberon   : mean-field theory included version           (2015) **/
/**       3.5 Miranda  : coupled-channels enhanced version            (2015) **/
/**                                                                          **/
/******************************************************************************/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <unistd.h>

using namespace std;
#define COH_TOPLEVEL

#include "structur.h"
#include "coh.h"
#include "global.h"
#include "parameter.h"
#include "output.h"
#include "coupling.h"
#include "elements.h"
#include "macs.h"
#include "terminate.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


static void     cohOutputOptions    (unsigned int, unsigned int);
static void     cohHelp             (void);
static void     cohAllocateMemory   (void);
static void     cohDeleteAllocated  (void);

/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/

//---------------------------------------
//      Defined in global.h

GlobalControl     ctl;               // booleans to control calculation
GlobalOption      opt;               // booleans to activate optional features
GlobalPrint       prn;               // booleans to control output
GlobalPrintExt    pex;               // booleans to control non-standard output
GlobalPrintMFT    pmf;               // booleabs to control mean-field theory output

//---------------------------------------
//      Extern variables defined in nucleus.h

CrossSection      crx;               // cross section data
GammaProduction   gml;               // gamma-ray line data 
Nucleus          *ncl;               // nucleus data
NuclearStructure *nst;               // nuclear structure data
Fission          *fbr;               // fission barrier data

//---------------------------------------
//      Extern variables defined in parameter.h

Adjustable       adj;                // adjustable parameters

//---------------------------------------
//      Defined in terminate.h

ostringstream     message;           // error or warning message string


/**********************************************************/
/*      Extern Variables                                  */
/**********************************************************/

#ifdef HAVE_CONFIG_H
string           version = VERSION;  // code version
#else
string           version = "ver.3.5";
#endif
double          *fact = NULL;        // factorials for coupling constants


/**********************************************************/
/*      Local Scope Variables inside coh.cpp              */
/**********************************************************/
static int       allocated_ncl = 0;  // total number of nucleus
static int       allocated_nst = 0;  // total number of structure data


/**********************************************************/
/*      CoH Starts Here                                   */
/**********************************************************/
int main(int argc, char *argv[])
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

//---------------------------------------
//      Command Line Options

  /*** default output options */
  unsigned int  propt = PRN_SYSTEM | PRN_XSECTION;
  unsigned int  prext = 0;

  /*** Lab energy, Z, and A numbers from command line */
  double        labE  = 0.0;
  string        elem  = "";
  int           targZ = 0;
  int           targA = 0;
  unsigned long nsim  = 1000;

  int           x;
  while((x=getopt(argc,argv,"p:q:e:a:z:m:fxbh"))!=-1){
    switch(x){
    case 'p': propt      = atoi(optarg);                  break;
    case 'q': prext      = atoi(optarg);                  break;
    case 'e': labE       = atof(optarg);                  break;
    case 'a': targA      = atoi(optarg);                  break;
    case 'z': elem       = optarg;
              targZ = element_getZ(&elem[0]);
              if(targZ == 0) {
                cerr << "ERROR     :unknown target name " << elem << endl;
                exit(-1);
              }                                           break;
    case 'm': ctl.montecarlo = ctl.exclusive  = true;
              nsim       = atoi(optarg);                  break;
    case 'f': ctl.fileout    = true;                      break;
    case 'x': ctl.exclusive  = true;                      break;
    case 'h': cohHelp();                                  break;
    case 'b': outBanner(); exit(0);                       break;
    case ':': cerr << "ERROR     :need a value for option " << x << endl;
              cohHelp();                                  break;
    case '?': cerr << "ERROR     :unknown option " << argv[argc-1] << endl;
              cohHelp();                                  break;
    }
  }


//---------------------------------------
//     Calc Flow Initialize

  /*** super-lazy mode */
  if(targZ > 0 && targA > 0) ctl.superlazy = true;


//---------------------------------------
//     Output Items

  /*** for exclusive spectrum calculation, all standard output will be turned off */
  if(ctl.exclusive){
    propt = 0;
    prext = prext & (PEX_PTABLE | PEX_DDX);
  }

  /*** set print-out global variable */
  cohOutputOptions(propt,prext);


//---------------------------------------
//     Call Main Loop

  /*** memory allocation */
  cohAllocateMemory();

  /*** main part */
  cohMainLoop(labE,targZ,targA,nsim);

  /***  delete allocated memory */
  cohDeleteAllocated();

  return 0;
}


/**********************************************************/
/*      Setting Output Options                            */
/**********************************************************/
void  cohOutputOptions(unsigned int p, unsigned int q)
{
  /*** if non-standard output option given, turn off standard outputs */
  if(q != 0) p = 0;

  /*** check and activate each bit */
  prn.system       = p & PRN_SYSTEM;
  prn.xsection     = p & PRN_XSECTION;
  prn.spectra      = p & PRN_SPECTRA;
  prn.levelexcite  = p & PRN_LEVELEXCITE;
  prn.gcascade     = p & PRN_GCASCADE;
  prn.population   = p & PRN_POPULATION;
  prn.transmission = p & PRN_TRANSMISSION;
  prn.angdist      = p & PRN_ANGDIST;
  prn.preeqparm    = p & PRN_PREEQPARM;
  prn.smatrix      = p & PRN_SMATRIX;

  pex.gbranch      = q & PEX_GBRANCH;
  pex.cumulevel    = q & PEX_CUMULEVEL;
  pex.density      = q & PEX_DENSITY;
  pex.fisdensity   = q & PEX_FISDENSITY;
  pex.ddx          = q & PEX_DDX;
  pex.ptable       = q & PEX_PTABLE;
  pex.gammaline    = q & PEX_GAMMALINE;
  pex.decaywidth   = q & PEX_DECAYWIDTH;
}


/**********************************************************/
/*      Help ! I need somebody                            */
/**********************************************************/
void cohHelp()
{
  outBanner();

  cerr << "coh : the Hauser Feshbach code  ver." << VERSION << endl;
  cerr << "coh [options] < input_file " << endl;
  cerr << "    -h          : this help" << endl;
  cerr << "    -e energy   :  incident energy override" << endl;
  cerr << "    -p N        :  N is sum of these numbers" << endl;
  cerr << "                   1 : system parameters" << endl;
  cerr << "                   2 : cross section" << endl;
  cerr << "                   4 : particle emission energy spectra" << endl;
  cerr << "                   8 : level excitation cross sections (inelastic scattering)" << endl;
  cerr << "                  16 : gamma-ray cascade" << endl;
  cerr << "                  32 : continuum bin population" << endl;
  cerr << "                  64 : transmission coefficients for the entrance channel" << endl;
  cerr << "                 128 : scattering angular distributions" << endl;
  cerr << "                 256 : preequilibrium parameters" << endl;
  cerr << "                 512 : S-matrix elements for the entrance channel" << endl;
  cerr << "    -q M        :  M is sum of these numbers" << endl;
  cerr << "                   1 : gamma-ray branching ratios" << endl;
  cerr << "                   2 : cumulative discrete levels" << endl;
  cerr << "                   4 : level densities" << endl;
  cerr << "                   8 : fission level densities" << endl;
  cerr << "                  16 : DDX data for ENDF formatting" << endl;
  cerr << "                  32 : compound nucleus decay probability matrix" << endl;
  cerr << "                  64 : all gamma lines including primary gamma" << endl;
  cerr << "    -f          : generate tabulated cross section output" << endl;
  cerr << "  ... or" << endl;
  cerr << "coh -e Ein -a targA -z targZ" << endl;
  cerr << "    for a superlazy mode for neutron-incident reactions" << endl;

  exit(0);
}


/**********************************************************/
/*      Allocate Memory                                   */
/**********************************************************/
void cohAllocateMemory()
{
  /*** Instantiation of total number of
       Nucleus class, DiscreteLevel, and Fission objects.
       For the Nucleus and NuclearStructure, only the first object will
       have the actual memory allocation in heap, and all the other
       objects will be allocated later dynamically. */
  try{
    ncl = new Nucleus [MAX_COMPOUND];
    ncl[allocated_ncl++].memalloc(MAX_LEVELS,MAX_ENERGY_BIN);

    /*** nuclear structure data */
    nst = new NuclearStructure [MAX_NUCLIDE];
    nst[allocated_nst++].memalloc(MAX_LEVELS);

    /*** fission barrier data */
    fbr = new Fission [MAX_FISS_CHANCE];
    for(int i=0 ; i<MAX_FISS_CHANCE ; i++){
      fbr[i].memalloc();
    }

    /*** cross section data */
    crx.memalloc();
      
    /*** gamma-ray line array */
    if(pex.gammaline) gml.memalloc(MAX_GAMMALINE);

    /*** factorial */
    factorial_allocate();

    /*** parameter adjustment */
    adj.memalloc(MAX_ADJUSTABLES);
  }
  catch(bad_alloc &e){
    cerr << "ERROR     :memory allocation error" << endl;
    exit(-1);
  }
}


/**********************************************************/
/*      Heap Memory Allocation for Nucleus Object         */
/**********************************************************/
void cohAllocateNucleus()
{
  /*** Each call of this subroutine increments the total number of
       Nucleus object, and allocate the heap memory */
  try{
    ncl[allocated_ncl++].memalloc(MAX_LEVELS,MAX_ENERGY_BIN);
  }
  catch(bad_alloc &e){
    message << "memory allocation error for nucleus " << allocated_ncl;
    cohTerminateCode("cohAllocateNucleus");
  }
  if(allocated_ncl == MAX_COMPOUND){
    message << "number of compounds exceed the maximum " << allocated_ncl;
    cohTerminateCode("cohAllocateNucleus");
  }
}


/**********************************************************/
/*      Heap Memory Allocation for DiscreteLevel Object   */
/**********************************************************/
void cohAllocateNuclearStructure()
{
  /*** Each call of this subroutine increments the total number of
       NuclearStructure object, and allocate memory in the heap */
  try{
    nst[allocated_nst++].memalloc(MAX_LEVELS);
  }
  catch(bad_alloc &e){
    message << "memory allocation error for nuclear structure data " << allocated_nst;
    cohTerminateCode("cohAllocateNuclearStructure");
  }
  if(allocated_nst == MAX_NUCLIDE){
    message << "number of unique nuclides exceed the maximum " << allocated_nst;
    cohTerminateCode("cohAllocateNuclearStructure");
  }
}


/**********************************************************/
/*      Reset Object Pointers                             */
/**********************************************************/
void cohResetAllocationPointer()
{
  /*** Before starting a new calculation,
       move the pointer to the first object. */
  allocated_nst = 1;
  allocated_ncl = 1;
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void cohDeleteAllocated()
{
  delete [] ncl;
  delete [] nst;
  delete [] fbr;

  if(pex.gammaline) gml.memfree();

  factorial_delete();
}


/**********************************************************/
/*     Emergency Stop                                     */
/**********************************************************/
int cohTerminateCode(string module)
{
  /*** Release global storage */
  cohDeleteAllocated();

  /*** Exit code */
  if(module == "") cerr << "ERROR     :" << message.str() << endl;
  else             cerr << "ERROR     : [" << module << "] " << message.str() << endl;
  exit(-1);
}


void cohNotice(string module)
{
  if(module == "NOTE"){
    cerr << "NOTICE    : " << message.str() << endl;
  }
  else{
    cerr << "WARNING   : [" << module << "] " << message.str() << endl;
  }
  message.str("");
}
