/******************************************************************************/
/*  omoutput.cpp                                                              */
/*        data output for optical model calculations                          */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "optical.h"
#include "omoutput.h"
#include "outformat.h"


/**********************************************************/
/*     Output optical potential parameters                */
/**********************************************************/
void outOMP(int ctl, Optical *omp)
{
  outSectionHead("OPTICAL MODEL POTENTIAL");

  if(ctl >= 0){
    cout << "# OMP" << setw(3) << ctl << "   ";
    cout << "   V         Vs        Wv        Ws        Vso       Wso";
  }
  else  
    cout << cline << "   V         Vs        Wv        Ws        Vso       Wso";

  if(omp->rc != 0.0) cout << "       Coul" << endl;
  else cout << endl;

  cout << "      r[fm]";
  outVal(10,5,omp->r0);   outVal(10,5,omp->r0s);
  outVal(10,5,omp->rv);   outVal(10,5,omp->rs);
  outVal(10,5,omp->rvso); outVal(10,5,omp->rwso); if(omp->rc != 0.0) outVal(10,2,omp->rc);
  cout << endl;

  cout << "      a[fm]";
  outVal(10,5,omp->a0);   outVal(10,5,omp->a0s);
  outVal(10,5,omp->av);   outVal(10,5,omp->as);
  outVal(10,5,omp->avso); outVal(10,5,omp->awso);
  cout << endl;

  if(ctl < 0){
    cout << "  Depth(E0)";
    outVal(10,5,omp->v1);   outVal(10,5,omp->vs1);
    outVal(10,5,omp->wv1);  outVal(10,5,omp->ws1);
    outVal(10,5,omp->vso1); outVal(10,5,omp->wso1);
    cout << endl;

    cout << "  Depth(E1)";
    outVal(10,5,omp->v2);   outVal(10,5,omp->vs2);
    outVal(10,5,omp->wv2);  outVal(10,5,omp->ws2);
    outVal(10,5,omp->vso2); outVal(10,5,omp->wso2);
    cout << endl;

    cout << "  Depth(E2)";
    outVal(10,5,omp->v3);   outVal(10,5,omp->vs3);
    outVal(10,5,omp->wv3);  outVal(10,5,omp->ws3);
    outVal(10,5,omp->vso3); outVal(10,5,omp->wso3);
    cout << endl;
  }
  cout << " Depth[MeV]";

  outVal(10,5,omp->volume.real());     outVal(10,5,omp->surface.real());
  outVal(10,5,omp->volume.imag());     outVal(10,5,omp->surface.imag());
  outVal(10,5,omp->spin_orbit.real()); outVal(10,5,omp->spin_orbit.imag());
  cout << endl;
}


void outOMPtable(int ctl, Optical *omp)
{
  if(ctl == 0) outVal(8,3,outRetrieveLabE());
  else cout << "        ";

  outVal(8,3,omp->volume.real());     outVal(8,3,omp->r0);   outVal(8,3,omp->a0);
  outVal(8,3,omp->surface.real());    outVal(8,3,omp->r0s);  outVal(8,3,omp->a0s);
  outVal(8,3,omp->volume.imag());     outVal(8,3,omp->rv);   outVal(8,3,omp->av);
  outVal(8,3,omp->surface.imag());    outVal(8,3,omp->rs);   outVal(8,3,omp->as);
  outVal(8,3,omp->spin_orbit.real()); outVal(8,3,omp->rvso); outVal(8,3,omp->avso);
  outVal(8,3,omp->spin_orbit.imag()); outVal(8,3,omp->rwso); outVal(8,3,omp->awso);
  outVal(8,2,omp->rc);
  cout << endl;
}

/**********************************************************/
/*     Output collective states                           */
/**********************************************************/
void outCoupledState(NuclearModel m, int n, LevelData *lev)
{
  if(m == vibration) outSectionHead("VIBRATIONAL STATE DATA");
  else               outSectionHead("ROTATIONAL STATE DATA");
  cout << cline << "   Energy       Spin   WaveNum Eout[MeV]   Coulomb   Sig0   ";
  if(m == vibration) cout << " Phonon";
  cout << endl;
  for(int i=0 ; i<n ; i++){
    cout << "      Level";
    outVal(10,5,lev[i].excitation);
    outVal(10,1,lev[i].spin);
    outVal(10,5,lev[i].wave_number);
    outVal(10,5,lev[i].energy);
    outVal(10,5,lev[i].coulomb);
    outVal(10,5,lev[i].coulomb_scat0);
    if(m == vibration){ cout << setw(5) << lev[i].phonon; }
    cout << endl;
  }
}


/**********************************************************/
/*     Output deformation parameters                      */
/**********************************************************/
void outDeformation(NuclearModel m, int n, double *beta)
{
  outSectionHead("DEFORMATION PARAMETERS");
  cout << cline << "    Lambda     Beta" << endl;
  if(m == rotation){
    for(int i=0 ; i<n/2 ; i++){
      cout << "       Beta";
      outVal(10,(i+1)*2); outVal(10,5,beta[i]);
      cout << endl;
    }
  }
  else{
    for(int i=0 ; i<n ; i++){
      cout << "       Beta";
      outVal(10,i+2); outVal(10,5,beta[i]);
      cout << endl;
    }
  }
}


/**********************************************************/
/*      Discrete Level Population by Direct               */
/**********************************************************/
void outLevelExcite(int n0, int n1, double ecms, double ex, LevelData *lev, double *cdir)
{
  outSectionHead("DISCRETE DIRECT CROSS SECTION");
  cout << cline << "    Ex[MeV]     Spin    Eout[MeV]  Sigma[mb]" << endl;

  for(int i=n0 ; i<n1 ; i++){
    if(cdir[i]<0.0) continue;
    cout << "  Direct" << setw(3) << i;
    outVal(11,5,lev[i].excitation);
    outVal( 9,1,lev[i].spin);
    cout << "  ";
    outVal(ecms-lev[i].excitation+ex);
    outVal(cdir[i]);
    cout << endl;
  }
}


/**********************************************************/
/*      Output R-Matrix and Strength Functions            */
/**********************************************************/
void outRmatrix(int n, complex<double> *s, complex<double> *r)
{
  outSectionHead("R MATRIX AND STRENGTH FUNCTION");

  cout << cline << "   Rr         Ri         R'       PoleStr.   Str. Func."<<endl;
  for(int j=0 ; j<=n ; j++){
    cout << setw(11) << j;
    cout.setf(ios::scientific, ios::floatfield);
    cout << setprecision(3);
    cout << setw(11) << r[j].real() << setw(11) << r[j].imag();
    cout << setw(11) << s[j].real() << setw(11) << r[j].imag()/PI;
    cout << setw(11) << s[j].imag() << endl;
  }
}


/**********************************************************/
/*      Output S-Matrix Elements                          */
/**********************************************************/
void outSmatrix(int n, int spin2, complex<double> *s)
{
  outSectionHead("S MATRIX ELEMENTS");

  cout << "#        L ";
  switch(spin2){
  case 0 : cout << "  Re(L)"; break;
  case 1 : cout << "  Re(L+1/2)  Im(L+1/2)  Re(L-1/2)  Im(L-1/2)"; break;
  case 2 : cout << "  Re(L+1)    Im(L+1)    Re(L)      Im(L)      Re(L-1)    Im(L-1)"; break;
  default: break;
  }
  cout << endl;

  for(int l=0 ; l<=n ; l++){
    cout << setw(11) << l;
    cout.setf(ios::scientific, ios::floatfield);
    cout << setprecision(3);
    for(int j=0 ; j<=2 ; j++){
      if(j <= spin2){
        cout << setw(11) << s[3*l+j].real();
        cout << setw(11) << s[3*l+j].imag();
      }
    }
    cout << endl;
  }
}


/**********************************************************/
/*      Transmission Coefficients                         */
/**********************************************************/
void outTransmission(int n, int spin2, double f, double *tran)
{
  outSectionHead("TRANSMISSION COEFFICIENTS");

  cout << "#        L ";
  switch(spin2){
  case 0 : cout << "   T(L)                          "; break;
  case 1 : cout << "   T(L+1/2)   T(L-1/2)           "; break;
  case 2 : cout << "   T(L+1)     T(L)       T(L-1)  "; break;
  default: break;
  }
  cout << "Sig Partial" << endl;

  double t = 0.0;
  for(int l=0 ; l<=n ; l++){
    cout << setw(11) << l;
    for(int j=0 ; j<=2 ; j++){
      if(j<=spin2) outVal(tran[3*l+j]);
      else cout << blank;
    }

    double s = 0.0;
    if(spin2==0)      s =  (2*l+1)*tran[l*3  ];
    else if(spin2==1) s =    (l+1)*tran[l*3  ] +   l   *tran[l*3+1];
    else if(spin2==2) s = ((2*l+3)*tran[l*3  ] +(2*l+1)*tran[l*3+1] + (2*l-1)*tran[l*3+2])/3.0;
    outVal(s*f);
    t += s*f;
    cout << endl;

    if(s == 0.0) break;
  }

  cout << blank << blank << blank << "        Sum";
  outVal(t);
  cout << endl;
}


/**********************************************************/
/*      Output S-Matrix Elements for Coupled-Channels     */
/**********************************************************/
void outCCSmatrix(int nc, int j, int p, Collective *col, CCdata *cdt, complex<double> *s)
{
  outSectionHead("COUPLED-CHANNELS S MATRIX ELEMENTS");

  cout << "#   JPi  Nc";
  outVal(9,1,j/2.0);
  cout << setw(2) << ((p>0) ? '+' : '-') << setw(11) << nc << endl;
  cout << "#  n   I   l   j" << endl;

  for(int i=0 ; i<nc ; i++){

    cout.setf(ios::fixed, ios::floatfield);
    cout << setprecision(1);

    cout << setw(4) << cdt[i].level;
    cout << setw(4) << fabs(col->lev[cdt[i].level].spin);
    cout << setw(4) << cdt[i].chn.l;
    cout << setw(4) << cdt[i].chn.j2 * 0.5;

    cout.setf(ios::scientific, ios::floatfield);
    cout << setprecision(5);
    for(int j=0 ; j<nc ; j++){
      int ij = i*nc+j;
      cout << ' ' << setw(12) << s[ij].real() << ' ' << setw(12) << s[ij].imag();
    }
    cout << endl;
  }
}
