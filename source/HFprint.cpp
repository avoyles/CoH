/******************************************************************************/
/*  HFprint.cpp                                                               */
/*        output calculated results                                           */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "HF.h"
#include "outformat.h"

static void   HFPrintSectionHead (const char *);
static double PRNPointDensity (double, double, Basis *, double *);
static double PRNHermite  (int, double);
static double PRNLaguerre (int, int, double);


/**********************************************************/
/*      Write Header of Each Section                      */
/**********************************************************/
void HFPrintSectionHead(const char *sec)
{
  const string s = "[ Hartree Fock BCS Calculation ]";
  int n1 = s.length();
  int n2 = strlen(sec);

  string bar1 = "#";

  for(int i=0 ; i<DisplayWidth-n1-1 ; i++) bar1 += ".";
  cout << bar1 << s << endl;

  string bar2 = "#";
  for(int i=0 ; i<DisplayWidth-n2-1 ; i++) bar2 += " ";
  cout << bar2 << sec << endl;
}


/**********************************************************/
/*      Print Internal Values                             */
/**********************************************************/
void HFMonitor(int i, double d, HFEnergy *e, Moment *m)
{
  cout << "#----" << endl;
  cout << "#" << setw(4) << i;
  cout << "  zcm fm      Q20 fm^2    Q30 b3/2  ";
  cout << "  Q40 b2      QN          Zmin      " << endl;
  cout << "#    ";
  outVal(12,m->zcm);
  outVal(12,m->q20);
  outVal(12,m->q30 * 1.0e-3);
  outVal(12,m->q40 * 1.0e-4);
  outVal(12,m->qn);
  outVal(12,m->zmin);
  cout << endl;

  cout << "#    ";
  cout << "  Kinetic     Volume      Surface   ";
  cout << "  SpinOrbit   Coulomb     Total MeV   epsilon" << endl;
  cout << "#    ";
  outVal(12, e->kinetic);
  outVal(12, e->volume);
  outVal(12, e->surface);
  outVal(12, e->spinorbit);
  outVal(12, e->coulomb);
  outVal(12, e->total);
  outVal(12, d);
  cout << endl;
}


/**********************************************************/
/*      Print BCS Parameters                              */
/**********************************************************/
void HFPrintBCS(BCSParameter *b)
{
  HFPrintSectionHead("BCS PARAMETERS");

  cout << "#               Lambda     Delta      Gap" << endl;
  cout << " Neutron   ";
  outVal(11,4,b[0].chemical_potential);
  outVal(11,4,b[0].pairing_energy);
  outVal(11,4,b[0].average_gap);
  nl();
  cout << " Proton    ";
  outVal(11,4,b[1].chemical_potential);
  outVal(11,4,b[1].pairing_energy);
  outVal(11,4,b[1].average_gap);
  nl();
}


/**********************************************************/
/*      Print Calculated Energies                         */
/**********************************************************/
void HFPrintEnergy(HFEnergy *e)
{
  HFPrintSectionHead("TOTAL ENERGY");

  cout << " Nucleon Masses                  ";  outVal(e->nucleon);   nl();
  cout << " Kinetic Energy                  ";  outVal(e->kinetic);   nl();
  cout << " Volume Energy                   ";  outVal(e->volume);    nl();
  cout << " Surface Energy                  ";  outVal(e->surface);   nl();
  cout << " Spin-Orbit Energy               ";  outVal(e->spinorbit); nl();
  cout << " Coulomb Energy                  ";  outVal(e->coulomb);   nl();
  cout << "                                 -----------"; nl();
  cout << " TOTAL                           ";  outVal(e->total);  nl();
}


/**********************************************************/
/*      Print Calculated Moments                          */
/**********************************************************/
void HFPrintMoment(Moment *m)
{
  HFPrintSectionHead("NUCLEAR MOMENTS");

  cout << cline;
  cout << "  zcm[fm]     Q20[fm^2]   Q30[b3/2] ";
  cout << "  Q40[b2]     QN          Zmin[fm]" << endl;
  cout << blank;
  outVal(12,m->zcm);
  outVal(12,m->q20);
  outVal(12,m->q30 * 1.0e-3);
  outVal(12,m->q40 * 1.0e-4);
  outVal(12,m->qn);
  outVal(12,m->zmin);
  cout << endl;
}


/**********************************************************/
/*      Print Calculated Potential                        */
/**********************************************************/
void HFPrintPotential(Basis *basis, OneBodyPotential *pot, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  HFPrintSectionHead("POTENTIAL SHAPE");

  cout << "# R[fm]       Z[fm]     ";
  cout << "  V0(n)       EffMass(n)  Vso(n)      V0(p)       EffMass(p)  Vso(p)"<<endl;
  for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;
    for(int i=0 ; i<fl->nr ; i++){
      cout << setw(12) << sqrt(fl->x[i]/basis->bp2);
      cout << setw(12) << fh->x[j]/basis->bz;
      cout << setw(12) << pot[0].central[i][j];
      cout << setw(12) << pot[0].effmass[i][j];
      cout << setw(12) << pot[0].spinorb[i][j];
      cout << setw(12) << pot[1].central[i][j];
      cout << setw(12) << pot[1].effmass[i][j];
      cout << setw(12) << pot[1].spinorb[i][j] << endl;
    }
    cout << endl;
  }
}


/**********************************************************/
/*      Print Single-Particle States                      */
/**********************************************************/
void HFPrintSPState(Basis *b, SPEnergy *s)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  HFPrintSectionHead("SINGLE PARTICLE SPECTRA");

  const double ecut = 10.0;
  const double vcut = 1.0e-10;

  cout << "#     neutron shell                 proton shell" << endl;
  cout << "#     Omg2 Energy[MeV]  v2          Omg2 Energy[MeV]  v2" << endl;
  for(int k=0 ; k<b->nv ; k++){
    cout << setw(4) << k;
    cout << setw(6) << b->omega[s[0].block[k]];
    cout << setw(12) << s[0].energy[k];
    cout << setw(12) << s[0].v2[k];
    cout << setw(6) << b->omega[s[1].block[k]];
    cout << setw(12) << s[1].energy[k];
    cout << setw(12) << s[1].v2[k] << endl;

    if( (s[0].energy[k] > ecut) && (s[1].energy[k] > ecut) ) break;
    if( (s[0].v2[k] < vcut) && (s[1].v2[k] < vcut) ) break;
  }
}


/**********************************************************/
/*      Print Nuclear Density                             */
/**********************************************************/
void HFPrintDensityCalc(int nr, int nz, double dr, double dz, Basis *basis, double *matV, GridData *rho)
{
  for(int i=0 ; i<nr ; i++){
    cerr << setw(6) << i*dr;
    for(int j=-nz/2 ; j<nz/2 ; j++){
      if(j == 0) continue;
      rho->p[i][j] = PRNPointDensity(dr * i, dz * j, basis, matV);
    }
  }
  cerr << endl;
}


void HFPrintDensity(int nr, int nz, double dr, double dz, GridData *rhon, GridData *rhop)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  HFPrintSectionHead("DENSITY");

  /*** at Z=0, take average of both adjucent mesh */
  for(int i=0 ; i<nr ; i++){
    rhon->p[i][0] = 0.5*(rhon->p[i][-1] + rhon->p[i][1]);
    rhop->p[i][0] = 0.5*(rhop->p[i][-1] + rhop->p[i][1]);
  }

  for(int i=-nr+1 ; i<nr ; i++){
    int k = abs(i);
    for(int j=-nz/2 ; j<nz/2 ; j++){
      cout << setw(12) << i*dr << setw(12) << j*dz;
      cout << setw(12) << rhon->p[k][j];
      cout << setw(12) << rhop->p[k][j];
      cout << setw(12) << rhon->p[k][j] + rhop->p[k][j] << endl;
    }
    cout << endl;
  }
}


/**********************************************************/
/*      rho at Given (R,Z)                                */
/**********************************************************/
double PRNPointDensity(double r, double z, Basis *basis, double *matV)
{
  const double rhocut = 1.0e-07;

  double zeta = basis->bz * z;
  double eta  = basis->bp2 * r*r;

  double sum = 0.0;

  for(int kb=0 ; kb<basis->nb ; kb++){
    for(int k1=basis->index[kb] ; k1<basis->index[kb]+basis->ndata[kb] ; k1++){

      int nr1 = basis->state[k1].getNr();
      int nz1 = basis->state[k1].getNz();
      int l1  = basis->state[k1].getL();
      int s1  = basis->state[k1].getS();

      for(int k2=basis->index[kb] ; k2<basis->index[kb]+basis->ndata[kb] ; k2++){

        int nr2 = basis->state[k2].getNr();
        int nz2 = basis->state[k2].getNz();
        int l2  = basis->state[k2].getL();
        int s2  = basis->state[k2].getS();

        if(s1 != s2) continue;
        if(l1 != l2) continue;

        int k0 = (k2<=k1) ? k1*(k1+1)/2 + k2 : k2*(k2+1)/2 + k1;
        double rho12 = matV[k0];
        if( abs(rho12)<rhocut ) continue;

        double yh = PRNHermite(nz1,zeta) * PRNHermite(nz2,zeta)
          / sqrt(PI * pow(2.0, (double)(nz1+nz2)) * exp(fact[nz1] + fact[nz2]));
        double yl = PRNLaguerre(nr1,l1,eta) * PRNLaguerre(nr2,l2,eta)
          * sqrt(exp(fact[nr1] + fact[nr2] - fact[nr1 + l1] - fact[nr2 + l2]));
 
        sum += rho12 * yh * yl * pow(eta, (double)l1);
      }
    }
  }

  double d = sum * exp(-zeta*zeta - eta) * basis->bz * basis->bp2/PI;

  return(d);
}


double PRNHermite(int n, double x)
{
  double fh = 0.0;

  if(n == 0) fh = 1.0;
  else{
    int p = n%2;

    if((x == 0.0) && (p == 0)){
      fh = ((n/2)%2 == 0 ? 1.0 : -1.0) * exp( fact[n] - fact[n/2] );
    }
    else{
      double s = 0.0;
      for(int m=0 ; m<=(n-p)/2 ; m++){
        s += (m%2 == 0 ? 1.0 : -1.0) * pow(2.0*x,(n-2.0*m)) / exp(fact[m] + fact[n-2*m]);
      }
      fh = s * exp(fact[n]);
    }
  }
  return(fh);
}


double PRNLaguerre(int n, int k, double x)
{
  double fl = 0.0;

  if(n == 0) fl = 1.0;
  else{
    double s = 0.0;
    if(x == 0.0){
      s = 1.0 / exp(fact[n] + fact[k]);
    }
    else{
      for(int m=0 ; m<=n ; m++){
        s += pow(-x,(double)m) / exp(fact[m] + fact[n-m] + fact[k+m]);
      }
    }
    fl = s * exp(fact[n+k]);
  }

  return(fl);
}
