/******************************************************************************/
/*  FRDMprint.cpp                                                             */
/*        Output FRDM Calculation Results                                     */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

using namespace std;

#include "FRDM.h"
#include "outformat.h"

static void FRDMPrintSectionHead(const char *);

/**********************************************************/
/*      Write Header of Each Section                      */
/**********************************************************/
void FRDMPrintSectionHead(const char *sec)
{
  const string s = "[ FRDM Calculation ]";
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
/*      FRDM System Parameters                            */
/**********************************************************/
void FRDMPrintShape(MFTSystem *sys){
  FRDMPrintSectionHead("SHAPE PARAMETERS");

  cout << cline << "   epsilon1   epsilon2   epsilon3";
  cout <<"   epsilon4   epsilon5   epsilon6      gamma" << endl;
  cout << " DeformParm";
  outVal(11,4,sys->getEps(1));
  outVal(11,4,sys->getEps(2));
  outVal(11,4,sys->getEps(3));
  outVal(11,4,sys->getEps(4));
  outVal(11,4,sys->getEps(5));
  outVal(11,4,sys->getEps(6));
  outVal(11,4,sys->getGamma()); nl();
}


/**********************************************************/
/*      Print Macroscopic Energy                          */
/**********************************************************/
void FRDMPrintMacroQuantity(double eps, double delta, double w0, double *b, double bs, double bv, double bw, double br)
{
  FRDMPrintSectionHead("MACRO-PART QUANTITY");

  cout << cline << " eps-bar    delta-bar  W0/W0(sph)" <<  endl;
  cout << blank;
  outVal(eps);
  outVal(delta);
  outVal(w0); nl();

  cout << cline << " B1         B2         B3         B4        ";
  cout << " Bs         Bv         Bw         Br" << endl;
  cout << blank;
  for(int i=0 ; i<4 ; i++) outVal(b[i]);
  outVal(bs);
  outVal(bv);
  outVal(bw);
  outVal(br); nl();
}


void FRDMPrintMacroEnergy(double *e)
{
  FRDMPrintSectionHead("MACROSCOPIC ENERGIES");

  int i = 0;
  cout << " Mass Exess                      ";  outVal(e[i++]); nl();
  cout << " Volume Energy                   ";  outVal(e[i++]); nl();
  cout << " Surface Energy                  ";  outVal(e[i++]); nl();
  //  cout << " Curvature Energy                ";  outVal(e[i++]); nl();
  //  cout << " A0 Energy                       ";  outVal(e[i++]); nl();
  i += 2;
  cout << " Coulomb Energy                  ";  outVal(e[i++]); nl();
  cout << " Volume Redistribution Energy    ";  outVal(e[i++]); nl();
  cout << " Coulomb Exchange Energy         ";  outVal(e[i++]); nl();
  cout << " Surface Redistribution Energy   ";  outVal(e[i++]); nl();
  cout << " Proton Formfactor Energy        ";  outVal(e[i++]); nl();
  cout << " Charge Asymmetry Energy         ";  outVal(e[i++]); nl();
  cout << " Wigner Term Energy              ";  outVal(e[i++]); nl();
  cout << " Pairing Energy                  ";  outVal(e[i++]); nl();
  cout << " Bound Electrons                 ";  outVal(e[i++]); nl();
}


/**********************************************************/
/*      Print Spherical Harmonics Expansion               */
/**********************************************************/
void FRDMPrintBeta(const int lmax, double **alm, double R0)
{
  const double eps = 1e-10;

  FRDMPrintSectionHead("SHAPE BETA EXPANSION");

  cout << "#   L   M    a(L,M)      beta(L)   " << endl;

  for(int l=0 ; l<=lmax ; l++){
    for(int m=-l ; m<=l ; m++){
      double beta = alm[l][m] / R0;
      if(abs(beta) > eps){
        cout << setw(5) << l << setw(4) << m;
        cout << " "; outVal(alm[l][m]);
        cout << " "; outVal(beta);  nl();
      }
    }
  }
}


/**********************************************************/
/*      Print Generated FRDM Potential                    */
/**********************************************************/
void FRDMPrintPotential(Basis *basis, OneBodyPotential *pot, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  FRDMPrintSectionHead("POTENTIAL SHAPE");

  cout << "# R[fm]       Z[fm]     ";
  cout << "  V0(n)       Vso(n)      V0(p)       Vso(p)"<<endl;
  for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;
    for(int i=0 ; i<fl->nr ; i++){
      cout << setw(12) << sqrt(fl->x[i]/basis->bp2);
      cout << setw(12) << fh->x[j]/basis->bz;
      cout << setw(12) << pot[0].central[i][j];
      cout << setw(12) << pot[0].spinorb[i][j];
      cout << setw(12) << pot[1].central[i][j];
      cout << setw(12) << pot[1].spinorb[i][j] << endl;
    }
    cout << endl;
  }
}


/**********************************************************/
/*      Print Single-Particle States                      */
/**********************************************************/
void FRDMPrintSPState(const bool bcsflag, Basis *b, SPEnergy *s, LNParameter *bcs)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);

  if(bcsflag) FRDMPrintSectionHead("SHIFTED SINGLE PARTICLE SPECTRA");
  else        FRDMPrintSectionHead("SINGLE PARTICLE SPECTRA");

  const double ecut = 10.0;
  const double vcut = 1.0e-10;

  cout << "#     neutron shell                     proton shell" << endl;
  cout << "#     Omg2    Energy[MeV]  v2          Omg2    Energy[MeV]  v2" << endl;
  for(int k=0 ; k<b->nv ; k++){

    double en = s[0].energy[k];
    double ez = s[1].energy[k];

    if(bcsflag){
      if( (bcs[0].n1 <= k) && (k<=bcs[0].n2) ){
        en += (4.0*bcs[0].lambda2 - bcs[0].strength)*s[0].v2[k];
      }
      if( (bcs[1].n1 <= k) && (k<=bcs[1].n2) ){
        ez += (4.0*bcs[1].lambda2 - bcs[1].strength)*s[1].v2[k];
      }
    }

    cout << setw(4) << k;
    cout << setw(6) << b->omega[s[0].block[k]];
    cout << setw(3) << (s[0].parity[k] > 0.0 ? " + " : " - ");
    cout << setw(12) << en;
    cout << setw(12) << s[0].v2[k];
    cout << setw(6) << b->omega[s[1].block[k]];
    cout << setw(3) << (s[1].parity[k] > 0.0 ? " + " : " - ");
    cout << setw(12) << ez;
    cout << setw(12) << s[1].v2[k] << endl;

    if( (s[0].energy[k] > ecut) && (s[1].energy[k] > ecut) ) break;
    if( (s[0].v2[k] < vcut) && (s[1].v2[k] < vcut) ) break;
  }
}


/**********************************************************/
/*      Print BCS Parameters                              */
/**********************************************************/
void FRDMPrintBCS(LNParameter *bcs)
{
  FRDMPrintSectionHead("BCS PARAMETERS");

  cout << "#             Strength     Lambda    Lambda2";
  cout << "      Delta    Delta-G  DeltaMacr" << endl;
  cout << " Neutron   ";
  outVal(11,4,bcs[0].strength);
  outVal(11,4,bcs[0].chemical_potential);
  outVal(11,4,bcs[0].lambda2);
  outVal(11,4,bcs[0].pairing_gap);
  outVal(11,4,bcs[0].pairing_eff);
  outVal(11,4,bcs[0].pairing_mac);
  nl();
  cout << " Proton    ";
  outVal(11,4,bcs[1].strength);
  outVal(11,4,bcs[1].chemical_potential);
  outVal(11,4,bcs[1].lambda2);
  outVal(11,4,bcs[1].pairing_gap);
  outVal(11,4,bcs[1].pairing_eff);
  outVal(11,4,bcs[1].pairing_mac);
  nl();
}


/**********************************************************/
/*      Print Shell Correction Parameters                 */
/**********************************************************/
void FRDMPrintSC(MFTSystem *sys, SCParameter *scp)
{
  FRDMPrintSectionHead("SHELL CORRECTION PARAMETERS");

  cout << "#               Lambda    Rho-bar sp-density   nucl./13" << endl;
  cout << " Neutron   ";
  outVal(11,4,scp[0].chemical_potential);
  outVal(11,4,scp[0].rho);
  outVal(11,4,scp[0].rho*2.0);
  outVal(11,4,sys->getN()/13.0);
  nl();
  cout << " Proton    ";
  outVal(11,4,scp[1].chemical_potential);
  outVal(11,4,scp[1].rho);
  outVal(11,4,scp[1].rho*2.0);
  outVal(11,4,sys->getZ()/13.0);
  nl();
}


/**********************************************************/
/*      Print Microscopic Energy                          */
/**********************************************************/
void FRDMPrintMicroEnergy(LNParameter *bcs, SCParameter *scp)
{
  FRDMPrintSectionHead("MICROSCOPIC ENERGIES");

  cout << "#              Pairing    Average    Ep-<Ep>" << endl;
  cout << " Neutron   ";
  outVal(11,4,bcs[0].energy_micro);
  outVal(11,4,bcs[0].energy_average);
  outVal(11,4,bcs[0].energy_micro - bcs[0].energy_average); nl();

  cout << " Proton    ";
  outVal(11,4,bcs[1].energy_micro);
  outVal(11,4,bcs[1].energy_average);
  outVal(11,4,bcs[1].energy_micro - bcs[1].energy_average); nl();
  nl();

  cout << "#             ShellCor    Average    Es-<Es>" << endl;
  cout << " Neutron   ";
  outVal(11,4,scp[0].energy_micro);
  outVal(11,4,scp[0].energy_average);
  outVal(11,4,scp[0].energy_micro - scp[0].energy_average); nl();

  cout << " Proton    ";
  outVal(11,4,scp[1].energy_micro);
  outVal(11,4,scp[1].energy_average);
  outVal(11,4,scp[1].energy_micro - scp[1].energy_average); nl();
}


/**********************************************************/
/*      Print Total Micro / Macro Energies                */
/**********************************************************/
void FRDMPrintTotalEnergy(FRDMEnergy *e)
{
  FRDMPrintSectionHead("TOTAL ENERGY");

  cout << "                                    deformed    spherical"; nl();
  cout << " macroscopic energy              ";  outVal(e->macro);
  cout << "  ";                                 outVal(e->spherical); nl();
  cout << " microscopic pair+shell energies ";  outVal(e->micro); nl();
  cout << " zero-point energy               ";  outVal(e->zero); nl();
  cout << "                                 -----------";  nl();
  cout << " TOTAL                           ";  outVal(e->total); nl();
}

