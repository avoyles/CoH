/******************************************************************************/
/*  FRDMfileout.cpp                                                           */
/*        Print FRDM results on a file                                        */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "FRDM.h"
#include "global.h"


/**********************************************************/
/*      Print FRDM Nuclear Shape                          */
/**********************************************************/
void FRDMWrite3DMesh(string filename, PolarCoordinate *p)
{
  bool polygon = false;
  ofstream fp;
  fp.open(&filename[0]);

  fp.setf(ios::scientific, ios::floatfield);
  fp << setprecision(4);

  if(polygon){
    for(int i = 0 ; i<p->nt-1 ; i++){
      for(int j = 0 ; j<p->np-1 ; j++){
        double x11 = p->getR(i,j) * p->getX(i,j);
        double y11 = p->getR(i,j) * p->getY(i,j);
        double z11 = p->getR(i,j) * p->getZ(i);

        double x12 = p->getR(i+1,j) * p->getX(i+1,j);
        double y12 = p->getR(i+1,j) * p->getY(i+1,j);
        double z12 = p->getR(i+1,j) * p->getZ(i+1);

        double x21 = p->getR(i,j+1) * p->getX(i,j+1);
        double y21 = p->getR(i,j+1) * p->getY(i,j+1);
        double z21 = p->getR(i,j+1) * p->getZ(i);

        double x22 = p->getR(i+1,j+1) * p->getX(i+1,j+1);
        double y22 = p->getR(i+1,j+1) * p->getY(i+1,j+1);
        double z22 = p->getR(i+1,j+1) * p->getZ(i+1);

        fp << setw(12) << x11 << setw(12) << y11 << setw(12) << z11;
        fp << setw(12) << x12 << setw(12) << y12 << setw(12) << z12;
        fp << setw(12) << x22 << setw(12) << y22 << setw(12) << z22;
        fp << setw(12) << x21 << setw(12) << y21 << setw(12) << z21 << endl;
      }
    }
  }
  else{
    for(int i = 0 ; i<p->nt ; i++){
      for(int j = 0 ; j<p->np ; j++){

        double x = p->getR(i,j) * p->getX(i,j);
        double y = p->getR(i,j) * p->getY(i,j);
        double z = p->getR(i,j) * p->getZ(i);

        fp << setw(12) << x;
        fp << setw(12) << y;
        fp << setw(12) << z << endl;
      }
      fp << endl;
    }
  }

  fp.close();
}


/**********************************************************/
/*      Print FRDM Microscopic Model Parameters           */
/**********************************************************/
void FRDMWriteParameters(const int z, const int a, LNParameter *bcs, SCParameter *scp)
{
  string fname = fileFRDMParameter;

  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp) return;

  fp.setf(ios::scientific, ios::floatfield);
  fp << setprecision(4);

  fp << setw(4) << z << setw(4) << a;
  for(int p=0 ; p<=1 ; p++){
    fp << setw(12) << bcs[p].pairing_gap;
    fp << setw(12) << bcs[p].lambda2;
    fp << setw(12) << bcs[p].chemical_potential;
    fp << setw(12) << scp[p].chemical_potential;
    fp << setw(12) << scp[p].rho * 2.0;
  }

  double ep[2],es[2];
  for(int p=0 ; p<=1 ; p++){
    ep[p] = bcs[p].energy_micro - bcs[p].energy_average;
    es[p] = scp[p].energy_micro - scp[p].energy_average;
  }

  fp << setw(12) << ep[0] << setw(12) << ep[1] << setw(12) << ep[0]+ep[1];
  fp << setw(12) << es[0] << setw(12) << es[1] << setw(12) << es[0]+es[1];

  fp << endl;
}


/**********************************************************/
/*      Print FRDM Energies                               */
/**********************************************************/
void FRDMWriteEnergy(string fname, MFTSystem *sys, FRDMEnergy *e)
{
  ofstream fp;
  fp.open(fname.c_str(),ios::out | ios::app);
  if(!fp) return;

  fp.setf(ios::scientific, ios::floatfield);
  fp << setprecision(4);

  fp << setw(4) << sys->getZ()
     << setw(4) << sys->getA();

  fp << setw(12) << sys->getEps(2)
     << setw(12) << sys->getEps(3)
     << setw(12) << sys->getEps(4);

  fp << setw(12) << e->spherical
     << setw(12) << e->macro
     << setw(12) << e->micro
     << setw(12) << e->total
     << endl;
}
