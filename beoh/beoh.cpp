/******************************************************************************/
/**                                                                          **/
/**   BeoH   :  Beta-Decay by CoH                                            **/
/**                                            Version 1.2 Naiad    (2017)   **/
/**                                                              T. Kawano   **/
/**   History                                                                **/
/**   1.0  2014 Apr. (Triton)    Imported from CoH3.3.2 Titania              **/
/**             Jul.             Combined with CoH3.4.0 Oberon               **/
/**   1.1  2017 Apr. (Nereid)    Simple Statistical Decay Mode               **/
/**   1.2  2017 Sep. (Naiad)     Prompt Fission Calculation Mode             **/
/**                                                                          **/
/******************************************************************************/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

using namespace std;
#define COH_TOPLEVEL

#include "beoh.h"
#include "beohoutput.h"
#include "global.h"
#include "parameter.h"
#include "coupling.h"
#include "elements.h"
#include "terminate.h"

int             main                (int, char *[]);
static void     cohOutputOptions    (unsigned int);
static void     cohHelp             (void);

static void     cohAllocateMemory   (void);
static void     cohDeleteAllocated  (void);

/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/

//---------------------------------------
//      Extern variables defined in nucleus.h

CrossSection      crx;                // cross section data
Nucleus          *ncl;                // nucleus data
NuclearStructure *nst;                // discrete level data
Fission          *fbr;                // fission barrier data

string            version = "ver.1.2 (2017 Sep)";


//---------------------------------------
//      Defined in global.h

GlobalControl     ctl;               // booleans to control calculation
GlobalOption      opt;               // booleans to activate optional features
GlobalPrint       prn;               // booleans to control output

//---------------------------------------
//      Extern variables defined in parameter.h

Adjustable        adj;               // adjustable parameters

//---------------------------------------
//      Defined in terminate.h

ostringstream     message;           // error or warning message string

//---------------------------------------
//      Extern variables defined in coupling.h

double           *fact;              // factorials for coupling constants


/**********************************************************/
/*      Local Scope Variables inside beoh.cpp             */
/**********************************************************/
static int        allocated_ncl = 0; // total number of nucleus
static int        allocated_nst = 0; // total number of nuclear structure data


/**********************************************************/
/*      BeoH Starts Here                                  */
/**********************************************************/
int main(int argc, char *argv[])
{

//---------------------------------------
//      Command Line Options

  /*** default output options */
  unsigned int  propt = PRN_SYSTEM | PRN_XSECTION | PRN_SPECTRA;


  /*** excitation energy, Z, and A numbers from command line */
  double        energy = 0.0, exwidth = 0.0, spinfact = 1.0, iniJ = 0.0, targE = 0.0, beta2 = 0.0;
  string        elem  = "";
  int           targZ = 0, targA = 0, iniP = 0;
  int           x;
  while((x = getopt(argc,argv,"p:e:d:a:z:j:J:f:k:b:h")) != -1){
    switch(x){
    case 'p': propt      = atoi(optarg);                  break;
    case 'e': energy     = atof(optarg);                  break;
    case 'd': exwidth    = atof(optarg);                  break;
    case 'a': targA      = atoi(optarg);                  break;
    case 'z': elem       = optarg;
              targZ = element_getZ(&elem[0]);
              if(targZ == 0) {
                cerr << "ERROR     :unknown target name " << elem << endl;
                exit(-1);
              }                                           break;
    case 'j': iniJ       = abs(atof(optarg));
              iniP       = (strstr(optarg,"-")) ? -1 : 1; break;
    case 'J': iniJ       = abs(atof(optarg));
              iniP       = -2;                            break;
    case 'f': spinfact   = atof(optarg);                  break;
    case 'k': targE      = atof(optarg);                  break;
    case 'b': beta2      = atof(optarg);                  break;
    case 'h': cohHelp();                                  break;
    case ':': cerr << "ERROR     :need a value for option " << x << endl;
              cohHelp();                                  break;
    case '?': cerr << "ERROR     :unknown option " << argv[argc-1] << endl;
              cohHelp();                                  break;
    }
  }

  /*** set print-out global variable */
  cohOutputOptions(propt);

  /*** memory allocation */
  cohAllocateMemory();

  /*** main part */
  if((targZ > 0) && (targA > 0) && (energy > 0.0)){
    beohCGMCompatibleMode(targZ, targA, iniJ, iniP, energy, exwidth, spinfact, targE, beta2);
  }
  else{
    beohMainLoop();
  }

  /***  delete allocated memory */
  cohDeleteAllocated();

  return 0;
}


/**********************************************************/
/*      Setting Output Options                            */
/**********************************************************/
void  cohOutputOptions(unsigned int p)
{
  /*** check and activate each bit */
  prn.system       = p & PRN_SYSTEM;
  prn.xsection     = p & PRN_XSECTION;
  prn.spectra      = p & PRN_SPECTRA;
}


/**********************************************************/
/*      Help                                              */
/**********************************************************/
void cohHelp()
{
  outBanner();

  cerr << endl;
  cerr << "beoh -e Excitation -a targA -z targZ [-p -j -f -b]" << endl;
  cerr << "    -p N        :  N is sum of these numbers" << endl;
  cerr << "                   1 : system parameters" << endl;
  cerr << "                   2 : ground state production probabilities" << endl;
  cerr << "                   4 : particle emission energy spectra" << endl;
  cerr << "    -j Spin_and_Parity" << endl;
  cerr << "                   spin J=3/2 -> 3.5 as input" << endl;
  cerr << "                   parity is given by its sign, +3.5, -0, etc" << endl;
  cerr << "    -f SpinFactor" << endl;
  cerr << "                   scaling factor to determine the initial spin distribution" << endl;
  cerr << "    -b Deformation_Parameter" << endl;
  cerr << "                   deformation parameter, beta2, if deformed" << endl;
  exit(0);
}


/**********************************************************/
/*      Allocate Memory                                   */
/**********************************************************/
void cohAllocateMemory()
{
  try{
    /*** compound nucleus, all other nuclei are allocated dynamically */
    ncl = new Nucleus [MAX_COMPOUND];
    ncl[allocated_ncl++].memalloc(MAX_LEVELS,MAX_ENERGY_BIN);

    /*** discrete level data */
    nst = new NuclearStructure [MAX_NUCLIDE];
    nst[allocated_nst++].memalloc(MAX_LEVELS);

    /*** fission barrier data */
    fbr = new Fission [MAX_FISS_CHANCE];
    for(int i=0 ; i<MAX_FISS_CHANCE ; i++){
      fbr[i].memalloc();
    }

    /*** residual nucleus production cross section */
    crx.prod = new ParticleCount [MAX_COMPOUND];

    /*** particle emission spectra */
    crx.spectalloc(MAX_CHANNEL+2,MAX_ENERGY_BIN);

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
/*      Instantiation Nucleus Object                      */
/**********************************************************/
void cohAllocateNucleus()
{
  try{
    ncl[allocated_ncl++].memalloc(MAX_LEVELS,MAX_ENERGY_BIN);
  }
  catch(bad_alloc &e){
    message << "memory allocation error for nucleus " << allocated_ncl;
    cohTerminateCode("cohAllocateNucleus");
  }
}


/**********************************************************/
/*      Instantiation DiscreteLevel Object                */
/**********************************************************/
void cohAllocateNuclearStructure()
{
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
/*      Reset Nucleus Object                              */
/**********************************************************/
void cohResetNucleus()
{
  delete [] ncl;
  allocated_ncl = 0;

  /*** re-allocate nucleus up to MAX_COMPOUND */
  ncl = new Nucleus [MAX_COMPOUND];
  ncl[allocated_ncl++].memalloc(MAX_LEVELS,MAX_ENERGY_BIN);
}


/**********************************************************/
/*      Reset DiscreteLevel Object                        */
/**********************************************************/
void cohResetNuclearStructure()
{
  delete [] nst;
  allocated_nst = 0;

  /*** re-allocate discrete level data up to MAX_NUCLIDE */
  nst = new NuclearStructure [MAX_NUCLIDE];
  nst[allocated_nst++].memalloc(MAX_LEVELS);
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void cohDeleteAllocated()
{
  delete [] ncl;
  delete [] nst;
  delete [] fbr;
  delete [] crx.prod;
  crx.spectfree();
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
