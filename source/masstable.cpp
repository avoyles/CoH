/******************************************************************************/
/*  masstable.cpp                                                             */
/*        find nuclear mass from Mass Table                                   */
/******************************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>

using namespace std;

#include "masstable.h"
#include "terminate.h"

//#include "masstable_aw95.h"             // Audi Wapstra 1995 table
//#include "masstable_ripl2.h"            // AW95 + FRDM95 from RIPL2
//#include "masstable_ripl3.h"            // AW03 + FRDM95 from RIPL3
//#include "masstable_audi2011.h"         // AW11 + FRDM95 from RIPL3
#include "masstable_audi2012_frdm2012.h"  // AW12 + FRDM2012

/*************************************************/
/*  Read Mass Data from File                     */
/*************************************************/
void massReadFile(string filename)
{
  ifstream fp;
  string   str;
  
  fp.open(&filename[0]);
  if(!fp){
    message << "file for mass data " << filename << " does not exist!";
    cohTerminateCode("massReadFile");
  }
  
  getline(fp,str);
  int n = atoi(str.c_str());
  int m = sizeof(MassTable)/sizeof(MassExcess);

  if(n > m) cerr << "too big mass data might be truncated" << endl;

  /*** overwrite existing mass table */
  unsigned int za;
  float mass;
  int i=0;
  while(!fp.eof()){
    fp >> za >> mass;
    MassTable[i].za   = za;
    MassTable[i].mass = mass;
    i++;

    if(i >= m){
      nMassTable = i;
      break;
    }
  }
  fp.close();
}  



/*************************************************/
/*  Mass Excess, Terminate Code If Not Found     */
/*************************************************/
double mass_excess(const int z, const int a)
{
  double    mx  = 0.0;
  unsigned int za = z*1000+a;

  bool found = false;
  for(int i=0 ; i<nMassTable ; i++){
    if(MassTable[i].za == za){
      found = true;
      mx = MassTable[i].mass;
      break;
    }
  }

  if(!found){
    message << "mass data for Z " << z << " - A " << a << " not found";
    cohTerminateCode("mass_exess");
  }

  return(mx);
}


/*************************************************/
/*  Mass Excess, Non-Stop Mode                   */
/*************************************************/
double mass_excess(const int z, const int a, bool *found)
{
  double    mx  = 0.0;
  unsigned int za = z*1000+a;

  *found = false;
  for(int i=0 ; i<nMassTable ; i++){
    if(MassTable[i].za == za){
      *found = true;
      mx = MassTable[i].mass;
      break;
    }
  }

  /*** huge value returned, so that separation energy will be infinity */
  if(! *found) mx = 1e+10;

  return(mx);
}

