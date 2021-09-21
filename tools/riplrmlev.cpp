/******************************************************************************/
/**                                                                          **/
/**   riplrmlev : Utility Program to Remove Levels from RIPL level Dababase  **/
/**                                                              T. Kawano   **/
/**                                                                          **/
/******************************************************************************/

#include <sstream>
#include <ostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

using namespace std;

#include "riplrmlev.h"

static inline double readDbl(string s, const int i0, const int i1)
{ string d = s.substr(i0,i1); return(atof(&d[0])); }

static inline int readInt(string s, const int i0, const int i1)
{ string d = s.substr(i0,i1); return(atoi(&d[0])); }

int main(int, char **);


/**********************************************************/
/*      Main Program                                      */
/**********************************************************/
int main(int argc, char *argv[])
{
  if(argc == 1){ cerr << "usage: riplrmlev z000.dat" << endl;  exit(-1); }


  /*** allocate discrete level data object */
  RIPLLevel RIPL;
  try{
    RIPL.memalloc(MAX_LEVELS,MAX_GAMMA_BRANCH);
  }
  catch(bad_alloc){ cerr << "memory allocation error" << endl;  exit(-1); }


  /*** open data file */
  ifstream fp;
  fp.open(argv[1]);
  if(!fp){ cerr << "data file cannot open" << endl; exit(-1); }


  /*** start reading the given file */
  string str;
  while(getline(fp,str)){

    /*** search for Z and A entry in the file */
    RIPL.name = str.substr(0,5);

    int a = readInt(str, 5,5);
    int z = readInt(str,10,5);
    RIPL.za.setZA(z,a);

    RIPL.nol  = readInt(str,15, 5);
    RIPL.nog  = readInt(str,20, 5);
    RIPL.nmax = readInt(str,25, 5);
    RIPL.nc   = readInt(str,30, 5);
    RIPL.sn   = readDbl(str,35,12);
    RIPL.sp   = readDbl(str,47,12);

    if(RIPL.nol >= MAX_LEVELS){
      cerr << "too many discrete levels for Z = " << z << "  A = " << a << endl;
      exit(-1);
    }

    /*** for all discrete levels */
    for(int i=0 ; i<RIPL.nol ; i++){
      getline(fp,str);
      RIPL.lev[i].energy   = readDbl(str, 4,10);
      RIPL.lev[i].spin     = readDbl(str,15, 5);
      RIPL.lev[i].parity   = readInt(str,20, 3);
      RIPL.lev[i].halflife = readDbl(str,25, 9);
      RIPL.lev[i].ngamma   = readInt(str,34, 3);
      RIPL.lev[i].info     = str.substr(37,95);

      RIPL.lev[i].flag = (str[0] == '#' ) ? 1 : 0;

      if(RIPL.lev[i].ngamma >= MAX_GAMMA_BRANCH){
        cerr << "too many gamma branch for Z = " << z << "  A = " << a << " at I = " << i << endl;
        exit(-1);
      }

      /*** for gamma-ray branches */
      for(int j=0 ; j<RIPL.lev[i].ngamma ; j++){
        getline(fp,str);
        RIPL.lev[i].fs[j] = readInt(str,39, 4) - 1;
        RIPL.lev[i].eg[j] = readDbl(str,43,11);
        RIPL.lev[i].pg[j] = readDbl(str,54,11);  // Pg = Pe/(1+ICC)
        RIPL.lev[i].pe[j] = readDbl(str,65,11);
        RIPL.lev[i].cc[j] = readDbl(str,76,11);
      }
    }

    /*** repeat until all removal flags disappear */
    while(RIPLCheckRemoveFlag(&RIPL)){

      /*** remove one of levels when #-mark is given */
      RIPLRemoveLevel(&RIPL);
    }

    /*** output current ZA data set */
    RIPLPrintLevel(&RIPL);
  }
  fp.close();

  return 0;
}


/**********************************************************/
/*      Main Process to Delete One Level                  */
/**********************************************************/
void RIPLRemoveLevel(RIPLLevel *r)
{
  /*** processing starts from the high side */
  int ix = 0;
  for(int i=r->nol-1 ; i>0 ; i--){
    if(r->lev[i].flag == 1){
      ix = i;
      break;
    }
  }

  /*** decrement number of gammas */
  if(r->lev[ix].ngamma > 0) r->nog -= r->lev[ix].ngamma;

  /*** adjust feeding from higher to removing level */
  for(int i=ix+1 ; i<r->nol ; i++){

    /*** check if there is a branch feeding to the removing level */
    int  jx = -1;
    for(int j=0 ; j<r->lev[i].ngamma ; j++){
      if(r->lev[i].fs[j] == ix){
        jx = j;
        break;
      }
    }
    if(jx == -1) continue;

    /*** save this branch */
    double pg = r->lev[i].pg[jx];
    double pe = r->lev[i].pe[jx];
    double cc = r->lev[i].cc[jx];

    /*** shift j-index */
    for(int j=jx+1 ; j<r->lev[i].ngamma ; j++){
      r->lev[i].fs[j-1] = r->lev[i].fs[j];
      r->lev[i].eg[j-1] = r->lev[i].eg[j];
      r->lev[i].pg[j-1] = r->lev[i].pg[j];
      r->lev[i].pe[j-1] = r->lev[i].pe[j];
      r->lev[i].cc[j-1] = r->lev[i].cc[j];
    }

    /*** append saved data */
    jx = r->lev[i].ngamma-1;
    r->lev[i].pg[jx] = pg;
    r->lev[i].pe[jx] = pe;
    r->lev[i].cc[jx] = cc;

    /*** check if there is a branch to the ground state */
    int jg = -1;
    for(int j=0 ; j<r->lev[i].ngamma ; j++){
      if(r->lev[i].fs[j] == 0){
        jg = j;
        break;
      }
    }

    /*** when no other gs transition exists, force decaying to gs */
    if(jg == -1){
      r->lev[i].fs[jx] = 0;
      r->lev[i].eg[jx] = r->lev[i].energy;
    }
    /*** combine transitions */
    else{
      r->lev[i].pg[jg] += r->lev[i].pg[jx];
      r->lev[i].pe[jg] += r->lev[i].pe[jx];
      r->lev[i].cc[jg]  = (r->lev[i].pg[jg] == 0.0) ? 0.0 : r->lev[i].pe[jg] / r->lev[i].pg[jg] - 1.0;
      r->lev[i].ngamma --;
    }
  }

  /*** first, save removing level */
  Level slev = r->lev[ix];

  /*** shift levels above ix */
  for(int i=ix+1 ; i<r->nol ; i++){
    r->lev[i-1] = r->lev[i];
  }

  /*** restore removed level at the top to keep memory pointers */
  r->lev[r->nol - 1] = slev;

  /*** decrement number of levels */
  r->nol --;
  if(r->nmax > ix) r->nmax --;
  if(r->nc > ix) r->nc --;

  /*** clearn the modifier flag */
  r->lev[ix].flag = 0;
}


/**********************************************************/
/*      Print Discrete Level Data for (Z,A)               */
/**********************************************************/
void RIPLPrintLevel(RIPLLevel *r)
{
  cout << setw(5) << r->name;
  cout << setw(5) << r->za.getA();
  cout << setw(5) << r->za.getZ();
  cout << setw(5) << r->nol;
  cout << setw(5) << r->nog;
  cout << setw(5) << r->nmax;
  cout << setw(5) << r->nc;

  cout.setf(ios::fixed, ios::floatfield);
  cout << setw(12) << setprecision(6) << r->sn;
  cout << setw(12) << setprecision(6) << r->sp;
  cout << endl;

  for(int i=0 ; i<r->nol ; i++){

    cout << setw(3) << i+1;

    cout.setf(ios::fixed, ios::floatfield);
    cout << setw(11) << setprecision(6) << r->lev[i].energy;
    cout << setw(6)  << setprecision(1) << r->lev[i].spin;
    cout << setw(3) << r->lev[i].parity;

    cout.setf(ios::scientific, ios::floatfield);

    if(r->lev[i].halflife == 0.0)
      cout << "           ";
    else
      cout << setw(11) << setprecision(2) << r->lev[i].halflife;

    cout << setw(3) << r->lev[i].ngamma;
    cout << r->lev[i].info;
    cout << endl;

    for(int j=0 ; j<r->lev[i].ngamma ; j++){
      cout << "                                       ";
      cout << setw(4) << r->lev[i].fs[j] + 1;
      cout << " ";

      cout.setf(ios::fixed, ios::floatfield);
      cout << setw(10) << setprecision(3) << r->lev[i].eg[j];

      cout.setf(ios::scientific, ios::floatfield);
      cout << setw(11) << setprecision(3) << r->lev[i].pg[j];
      cout << setw(11) << setprecision(3) << r->lev[i].pe[j];
      cout << setw(11) << setprecision(3) << r->lev[i].cc[j];
      cout << endl;
    }
  }
}


/**********************************************************/
/*      Check Whethe Removal Flag Is Given                */
/**********************************************************/
bool RIPLCheckRemoveFlag(RIPLLevel *r)
{
  bool f = false;
  for(int i=0 ; i<r->nol ; i++){
    if(r->lev[i].flag == 1){
      f = true;
      break;
    }
  }
  return f;
}

