/******************************************************************************/
/*  output.cpp                                                                */
/*        main data output                                                    */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

#include "physicalconstant.h"
#include "structur.h"
#include "output.h"
#include "parameter.h"
#include "nucleus.h"
#include "elements.h"
#include "global.h"


/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/

extern string version;

static string particle_name[9]={
   "      gamma","    neutron","     proton","      alpha",
   "   deuteron","     triton","     helion","    fission",
   "    unknown"};

static string p_name = "gnpadth";

static double labE = 0.0;   // laboratory energy
static double sigR = 0.0;   // total reaction cross section
static ZAnumber targZA(0,0);

void outSetSigmaReaction(const int n)
{
  sigR = 0.0;
  for(int i=0 ; i<n ; i++) if(crx.prod[i].xsec > 0.0) sigR += crx.prod[i].xsec;
}

#include "outbanner.h"
#include "outformat.h"


/**********************************************************/
/*      Header                                            */
/**********************************************************/
void outTitle(char *str)
{
#ifdef PRINT_BANNER
  outBanner1();
#endif
  outSectionHead(&version[0]);
  cout <<cline;
  cout << &str[10] << endl;
}


/**********************************************************/
/*      Print Main Banner                                 */
/**********************************************************/
void outBanner()
{
  static int val[] = {0,0,2,0,2,2,4,5,6,7,8,9,10};
  time_t c;

  /*** get today */
  time(&c);
  struct tm *p = localtime(&c);
  int year  = p->tm_year + 1900;
  int month = p->tm_mon  + 1;
  int day   = p->tm_mday;

  /*** calculate moon age */
  int a = ( (year - 11)%19 ) * 11;
  int b = val[month];
  int x = (a + b + day)%30;

  outSectionHead(&version[0]);
  if(x < 10)       cout <<banner1 << endl;
  else if(x < 20)  cout <<banner2 << endl;
  else             cout <<banner3 << endl;

  string bar = "#";
  for(int i=0 ; i<DisplayWidth-12 ; i++) bar += "#";
  cout << bar << " Moonage" << setw(3) << x << endl;
}


/**********************************************************/
/*      Write Z and A                                     */
/**********************************************************/
void outZA(ZAnumber *za)
{
  char element[3];

  if(za->getZ() >= N_ELEMENTS){
    element[0] = element[1] = ' ';
  }
  else{
    strncpy(element,element_name[za->getZ()].c_str(),2);
    if(strlen(element_name[za->getZ()].c_str()) == 1) element[1] = ' ';
  }
  element[2] = '\0';

  cout << "  ";
  cout << setw(3) << setfill('0') << za->getZ() << '-';
  cout << setw(3) << setfill('0') << za->getA();
  cout << setw(2) << element;
  cout << setfill(' ');
}


/**********************************************************/
/*      Write Header of Each Section                      */
/**********************************************************/
void outSectionHead(const char *sec)
{
  int n = strlen(sec);

  string bar1 = "#";

  if(labE == 0.0){
    for(int i=0 ; i<DisplayWidth-1 ; i++) bar1 += "#";
    cout << bar1 << endl;
  }else{
    for(int i=0 ; i<DisplayWidth-27 ; i++) bar1 += ".";
    cout << bar1 << "/";  outZA(&targZA);  cout << " // ";  outVal(10,6,labE); nl();
  }

  string bar2 = "#";
  for(int i=0 ; i<DisplayWidth-n-1 ; i++) bar2 += " ";
  cout << bar2 << sec << endl;
}


/**********************************************************/
/*      System Parameters                                 */
/**********************************************************/
void outSystem(System *sys, bool flag)
{
  /*** save Z, A, E for the first time */
  targZA.setZA(sys->target.getZ(),sys->target.getA());
  labE  = sys->lab_energy;
  if(!flag) return;

  outSectionHead("SYSTEM PARAMETERS");
  cout << cline << "  AtomicNum    MassNum" << endl;
  cout << "#    Target"; outVal(11,sys->target.getZ()); outVal(11,sys->target.getA()); nl();
  cout << "#  Incident"; outVal(11,sys->incident.za.getZ()); outVal(11,sys->incident.za.getA()); nl();
  cout << cline << "  Ecms[MeV]  Elab[MeV] ReduceMass WaveNumber" << endl;
  cout << "     Energy"; outVal(sys->cms_energy);   outVal(sys->lab_energy);
                         outVal(sys->reduced_mass); outVal(sys->wave_number); nl();
}


/**********************************************************/
/*      Target State                                      */
/**********************************************************/
void outTargetState(const int targid, const int targlev, double ex)
{
  outSectionHead("TARGET STATE");
  cout << cline << " Excitation  Spin" << endl;

  outZA(&ncl[targid].za);
  outVal(ex);
  outVal(5,1,ncl[targid].lev[targlev].spin);
  char p = (ncl[targid].lev[targlev].parity < 0) ? '-' : '+';
  cout << p << endl;
}


/**********************************************************/
/*      Compound Nucleus, Z and A, Binding Energies       */
/**********************************************************/
void outCompound(const int n, Pdata *pdt)
{
  outSectionHead("COMPOUND NUCLEUS DATA");
  cout << cline << "    Z     A  MassExess    Ex(max) Separation     N-cont    N-level" << endl;

  for(int i=0 ; i<n ; i++){
    outZA(&ncl[i].za);
    cout << blank;
    outVal(11,4,ncl[i].mass_excess); outVal(11,4,ncl[i].max_energy);
    cout << blank;
    outVal(11,ncl[i].ncont); outVal(11,ncl[i].ndisc); nl();
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(ncl[i].cdt[j].status){
        cout << particle_name[j];
        outZA(&ncl[ ncl[i].cdt[j].next ].za);
        outVal(11,4,pdt[j].mass_excess);
        cout << blank;
        outVal(11,4,ncl[i].cdt[j].binding_energy); nl();
      }
    }
  }
}


/**********************************************************/
/*      Level Density Parameters of Each Compound         */
/**********************************************************/
void outLevelDensity(const int n, double d0)
{
  outSectionHead("LEVEL DENSITY PARAMETERS");
  if(d0 > 0.0)
    cout << cline <<"       a       Pairing     T          E0         Em         Sigma2     Eshell         D0" << endl;
  else
    cout << cline <<"       a       Pairing     T          E0         Em         Sigma2     Eshell" << endl;
  for(int i=0 ; i<n ; i++){
    if(ncl[i].ncont > 0){
      /*** check redundancy */
      bool found = false;
      for(int j=0 ; j<i ; j++){
        if(ncl[i].za == ncl[j].za){
          found = true;
          break;
        }
      }
      if(found) continue;

      outZA(&ncl[i].za);
      outVal(11,4,ncl[i].ldp.a);
      outVal(11,4,ncl[i].ldp.pairing_energy);
      outVal(11,4,ncl[i].ldp.temperature);
      outVal(11,4,ncl[i].ldp.E0);
      outVal(11,4,ncl[i].ldp.match_energy);
      outVal(11,4,ncl[i].ldp.sigma0*ncl[i].ldp.sigma0);
      outVal(11,6,ncl[i].ldp.shell_correct);
      if( (i == 0) && (d0 > 0.0) ) outVal(d0);
      nl();
    }
  }
}


/**********************************************************/
/*      GDR Parameter                                     */
/**********************************************************/
void outGDR(GDR *gdr, double gg)
{
  outSectionHead("GIANT DIPOLE RESONANCE DATA");
  if(gg > 0.0){
    cout << cline << "Energy[MeV] Width[MeV] Sigma0[mb] <Gam>[MeV]" << endl;
  }
  else{
    cout << cline << "Energy[MeV] Width[MeV] Sigma0[mb]" << endl;
  }
  for(int i=0 ; i<MAX_GDR ; i++){
    if(gdr[i].getEnergy() > 0.0){
      cout << "         " << gdr[i].getEM() << setw(1) << gdr[i].getL();
      outVal(11,4,gdr[i].getEnergy());
      outVal(11,4,gdr[i].getWidth());
      outVal(11,6,gdr[i].getSigma()); nl();
    }
  }
  if(gg > 0.0){
    cout << "   Average Gamma Width" << blank << blank;
    outVal(gg); nl();
  }

}


/**********************************************************/
/*      Fission Barriers                                  */
/**********************************************************/
void outFissionBarrier(const int n)
{
  outSectionHead("FISSION BARRIER PARAMETERS");

  cout << cline << blank << "Height[MeV] Width[MeV] KBand[MeV]       Elmax[MeV]" << endl;
  for(int i=0 ; i<n ; i++){
    for(int m=0 ; m<MAX_HUMP ; m++){
      if(ncl[i].fissile){
        if(ncl[i].fission->barrier[m].height > 0.0){
          outZA(&ncl[i].za);
          if     (m == 0) cout << " 1stBarrier";
          else if(m == 1) cout << " 2ndBarrier";
          else if(m == 2) cout << " 3rdBarrier";
          else            cout << blank;
          outVal(11,4,ncl[i].fission->barrier[m].height);
          outVal(11,4,ncl[i].fission->barrier[m].curvature);
          cout << blank << "      ";
          outVal(ncl[i].fission->barrier[m].elmax); nl();

          for(int k=0 ; k<ncl[i].fission->barrier[m].nband ; k++){
            cout << blank << blank << blank << blank;
            outVal(ncl[i].fission->barrier[m].kband[k].excitation);
            outVal(5,1,ncl[i].fission->barrier[m].kband[k].k2/2.0);
            char p = (ncl[i].fission->barrier[m].kband[k].parity < 0) ? '-' : '+';
            cout << p << endl;

            /*** to print band levels */
#ifdef PRINT_BANDLEVELS
            int l0 = ((ncl[i].fission->barrier[m].kband[k].k2 == 0)
                   && (ncl[i].fission->barrier[m].kband[k].parity < 0)) ? 1 : 0;
            int ls = (ncl[i].fission->barrier[m].kband[k].k2 == 0) ? 2 : 1;
            double xk = ncl[i].fission->barrier[m].kband[k].k2/2.0;
            double hb = ncl[i].fission->barrier[m].inertia;
            double e0 = ncl[i].fission->barrier[m].kband[k].excitation;
            for(int l=l0 ; ; l+=ls){
              double xj = xk + l;
              double e = (xj*(xj+1.0)-xk*(xk+1.0))* hb + e0;
              if(e > ncl[i].fission->barrier[m].elmax) break;
              cout << blank << blank << blank << blank << blank;
              outVal(5,1,xj);
              cout << p;
              outVal(e); nl();
            }
#endif
          }
        }
      }
    }
  }
  if(ncl[0].fissile && ncl[0].fission->potwell[0].height > 0.0){
    nl();
    cout << cline << blank << "Height[MeV] Width[MeV] Aborb[MeV]" << endl;
    for(int i=0 ; i<n ; i++){
      for(int m=0 ; m<MAX_HUMP-1 ; m++){
        if(ncl[i].fissile){
          if(ncl[i].fission->potwell[m].height > 0.0){
            outZA(&ncl[i].za);
            if     (m == 0) cout << " Class-II  ";
            else if(m == 1) cout << " Class-III ";
            else            cout << blank;
            outVal(11,4,ncl[i].fission->potwell[m].height);
            outVal(11,4,ncl[i].fission->potwell[m].curvature);
            outVal(11,4,ncl[i].fission->potwell[m].absorb); nl();
          }
        }
      }
    }
  }
}


/**********************************************************/
/*      Cross Section                                     */
/**********************************************************/
void outCrossSection(const int ctl, double s1, double s2, double s3)
{
  outSectionHead("TOTAL REACTION CROSS SECTION");
  if(ctl == 0){
    cout << cline << "      Total    Elastic   Reaction" << endl;
    cout << "  Sigma[mb]";
    outVal(s1); outVal(s2); outVal(s3); nl();
  }
  else if(ctl == 1){
    cout << "  Sigma[mb]"; outVal(s3); nl();
  }
  else if(ctl == 2){
    cout << cline << "        GDR    Quasi-d   Reaction" << endl;
    cout << "  Sigma[mb]";
    outVal(s1); outVal(s2); outVal(s3); nl();
  }
}


/**********************************************************/
/*      Reaction Cross Section                            */
/**********************************************************/
void outReaction(const int n, const int targid, double ex)
{
  outSectionHead("INDIVIDUAL GROUND STATE PRODUCTION");
  cout << cline << "   Residual  Prod.[mb]  Qval[MeV]   Particles" << endl;
  cout << cline << blank << blank << blank << "  ";
  for(int j=1 ; j<(int)p_name.length() ; j++){
    cout << setw(2)<< p_name.substr(j,1);
  }
  nl();

  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].xsec > 0){
      cout << " Production";
      outZA(&ncl[i].za);
      outVal(lowfilter(crx.prod[i].xsec));
      outVal(11,5, lowfilter(ex + ncl[i].max_energy - ncl[targid].max_energy));
      cout << "  ";

      for(int j=1 ; j<MAX_CHANNEL ; j++) cout << setw(2) << (int)crx.prod[i].par[j];
      nl();
    }
  }
  cout << blank << "        sum";
  outVal(lowfilter(sigR)); nl();
}


/**********************************************************/
/*      Total Residual Nucleus Production Cross Section   */
/**********************************************************/
void outTotalResidual(CumulativeResidualProduct *res)
{
  outSectionHead("STABLE AND LONG-LIVED STATE PRODUCTIONS");
  if(res->getNcurrent() == 0) return;

  cout << "#       Residual    Ex[MeV]     Jpi     T(1/2)[s]  Prod.[mb] Meta" << endl;

  int c = 1;
  double sum = 0.0;
  for(unsigned int z=res->zmin ; z<=res->zmax ; z++){
    for(unsigned int a=res->amin ; a<=res->amax ; a++){
      ZAnumber za(z,a);
      for(int i=0 ; i<res->getNcurrent() ; i++){
        if(res->rp[i].za != za) continue;

        outVal(5,c++);
        outZA(&res->rp[i].za);
        cout << "  ";
        outVal(9,5,res->rp[i].energy);
        outVal(7,1,res->rp[i].spin);
        char p0 = (res->rp[i].parity < 0) ? '-' : '+';
        cout << p0 << "   ";
        if(res->rp[i].halflife < 0.0) cout << dashl;
        else outVal(res->rp[i].halflife);
        outVal(lowfilter(res->rp[i].production));
        if(res->rp[i].metaflag > 0) cout << setw(5) << res->rp[i].metaflag;
        nl();

        if(res->rp[i].metaflag == 0) sum += res->rp[i].production;
      }
    }
  }
  cout << cline << blank << blank << "         g.s.sum";
  outVal(lowfilter(sum)); nl();
  nl(); nl();
}


/**********************************************************/
/*      Isomeric Ratios for Long-Lived Nuclides           */
/**********************************************************/
void outIsomericRatio(CumulativeResidualProduct *res)
{
  outSectionHead("ISOMERIC RATIOS");
  if(res->getNcurrent() == 0) return;

  cout << "#       Residual    Ex[MeV]     Jpi     T(1/2)[s] Ratio      m/(g-m)" << endl;

  int c = 1;
  for(unsigned int z=res->zmin ; z<=res->zmax ; z++){
    for(unsigned int a=res->amin ; a<=res->amax ; a++){
      ZAnumber za(z,a);

      /*** look at isomer */
      for(int i1=0 ; i1<res->getNcurrent() ; i1++){
        if(res->rp[i1].za != za) continue;
        if(res->rp[i1].metaflag == 0) continue;

      /*** look for ground state of the same ZA */
        for(int i0=0 ; i0<res->getNcurrent() ; i0++){
          if((res->rp[i0].za == za) && (res->rp[i0].metaflag == 0)){

            /*** production of gs includes meta-production */
            double d  = res->rp[i0].production - res->rp[i1].production;
            double r0 = (res->rp[i0].production == 0.0) ? 0.0 : res->rp[i1].production / res->rp[i0].production;
            double r1 = (d == 0.0) ? 0.0 : res->rp[i1].production / d;

            outVal(5,c++);
            outZA(&za);
            cout << "  ";
            outVal(9,5,res->rp[i1].energy);
            outVal(7,1,res->rp[i1].spin);
            char p0 = (res->rp[i1].parity < 0) ? '-' : '+';
            cout << p0 << "   ";
            outVal(res->rp[i1].halflife);
            outVal(lowfilter(r0));
            outVal(lowfilter(r1));
            nl();
          }
        }
      }
    }
  }
  nl(); nl();
}


/**********************************************************/
/*      Fission Cross Section                             */
/**********************************************************/
void outFission(const int n)
{
  outSectionHead("FISSON CROSS SECTION");
  cout << cline << "   Nucleus Fission[mb] Fiss.Prob." << endl;

  double sum = 0.0;
  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].fiss > 0) sum += crx.prod[i].fiss;
  }

  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].fiss > 0){
      int c = crx.prod[i].par[neutron] + 1;
      cout << " FisChance" << setw(1) << c;
      outZA(&ncl[i].za);
      outVal(lowfilter(crx.prod[i].fiss));
      outVal(11,3,crx.prod[i].fiss/sum*100.0);
      nl();
    }
  }
  cout << blank << " TotFission";
  outVal(lowfilter(sum)); nl();
  cout << blank << "  Fiss+Prod";
  outVal(lowfilter(sum+sigR)); nl();
}


/**********************************************************/
/*      Particle Production Cross Section                 */
/**********************************************************/
void outParticleProduction(const int n, Channel *cdt, double **spc)
{
  double sum[MAX_CHANNEL], ave[MAX_CHANNEL];

  for(int j=0 ; j<MAX_CHANNEL ; j++) sum[j] = ave[j] = 0.0;

  /*** for particles, use the ground-state production */
  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].xsec > 0){
      for(int j=1 ; j<MAX_CHANNEL ; j++){
        sum[j] += (int)crx.prod[i].par[j] * crx.prod[i].xsec;
      }
    }
  }

  /*** for gamma-ray, use the emission spectrum */
  int kmax = (int)(ncl[0].max_energy/ncl[0].de) +1;
  for(int k=0 ; k<kmax ; k++){
    sum[0] += spc[0][k] * ncl[0].de;
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(cdt[j].status) ave[j] += spc[j][k] * ncl[0].de * k * ncl[0].de;
    }
  }

  
  outSectionHead("PARTICLE PRODUCTION");
  cout << cline << "   Particle  Prod.[mb] Multiplic. Emean[MeV]" << endl;

  cout << " ParticlPrd" << particle_name[0];
  outVal(lowfilter(sum[0]));
  cout << dashl;
  outVal(lowfilter(ave[0] / sum[0])); nl();

  for(int j=1 ; j<MAX_CHANNEL ; j++){
    if (!cdt[j].status) continue;

    /*** particle multiplicity and average energy */
    double mp = sum[j] / sigR;
    double em = (sum[j] == 0.0) ? 0.0 : ave[j] / sum[j];

    cout << " ParticlPrd" << particle_name[j];
    outVal(lowfilter(sum[j])); outVal(lowfilter(mp)); outVal(lowfilter(em)); nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Scattering Angular Distribution                   */
/**********************************************************/
void outAngularDistribution(const int ctl, int n0, int np, int step, ZAnumber *za)
{
  if(crx.theta[0] == 0.0) return;

  const int column = 8;
  int page = 1, m = 0;

  switch(ctl){
  case 0:
    outSectionHead("ELASTIC SCATTERING ANGULAR DISTRIBUTION");
    break;
  case 1:
    outSectionHead("CC ELASTIC AND INELASTIC SCATTERING ANGULAR DISTRIBUTION");
    break;
  case 2:
    outSectionHead("DWBA INELASTIC SCATTERING ANGULAR DISTRIBUTION");
    break;
  case 3:
    outSectionHead("COMPOUND PLUS DIRECT SCATTERING ANGULAR DISTRIBUTION");
    if(np > MAX_ANGDISTLEVELS) np=MAX_ANGDISTLEVELS;
    for(int k=np-1 ; k>=0 ; k--){
      if(crx.angdist[k][0] > 0.0){
        np = k+1;
        break;
      }
    }
    break;
  default:
    break;
  }
  
  if(np > 0) page = (np-n0-1)/column+1;

  for(int p=0 ; p<page ; p++){
    m = (np-n0>=1) ? min(np-n0-1,(p+1)*column-1) : 0;

    cout << cline; outZA(za);  cout << endl;
    cout << "#  Ang[deg]";
    for(int j=p*column+n0 ; j<=m+n0 ; j++) cout << setw(3) << j << " [mb/sr]";
    nl();

    for(int i=0 ; i<MAX_ANGDIST ; i++){
      int k = (int)crx.theta[i];
      if(k%step == 0){
        outVal(crx.theta[i]);
        for(int j=p*column+n0 ; j<=m+n0 ; j++) outVal(crx.angdist[j][i]);
        nl();
      }
    }
    nl();
    nl();
  }
}


/**********************************************************/
/*      Legendre Coefficients                             */
/**********************************************************/
void outLegendreCoefficient(const int ctl, int n0, int np, int id, ZAnumber *za)
{
  if(crx.theta[0] == 0.0) return;

  const int column = 8;
  int page = 1, m = 0;

  switch(ctl){
  case 0:
    outSectionHead("ELASTIC SCATTERING LEGENDRE COEFFICIENTS");
    break;
  case 1:
    outSectionHead("CC ELASTIC AND INELASTIC LEGENDRE COEFFICIENTS");
    break;
  case 2:
    outSectionHead("DWBA INELASTIC LEGENDRE COEFFICIENTS");
    break;
  case 3:
    outSectionHead("COMPOUND PLUS DIRECT REACTION LEGENDRE COEFFICIENTS");
    if(np > MAX_ANGDISTLEVELS) np=MAX_ANGDISTLEVELS;
    for(int k=np-1 ; k>=0 ; k--){
      if(crx.angdist[k][0] > 0.0){
        np = k+1;
        break;
      }
    }
    break;
  default:
    break;
  }
  
  if(np > 0) page = (np-n0-1)/column+1;

  for(int p=0 ; p<page ; p++){
    if(np-n0 >= 1) m = min(np-n0-1,(p+1)*column-1);

    cout << cline; outZA(za);  cout << endl;
    cout << "#        L ";
    for(int j=p*column+n0 ; j<=m+n0 ; j++) cout << setw(3) << j << "th level";
    nl();

    int jmax = 0;
    for(int j=MAX_J-1 ; j>=0 ; j--){
      bool flag = false;
      for(int k=p*column+n0 ; k<=m+n0 ; k++){
        if(crx.legcoef[id][k][j] != 0.0){ flag = true; break; }
      }
      if(flag){ jmax = j+1; break; }
    }

    cout.setf(ios::scientific, ios::floatfield);

    for(int j=0 ; j<=jmax ; j++){
      outVal(10,j);
      cout << " ";
      for(int k=p*column+n0 ; k<=m+n0 ; k++){
        double x0 = crx.legcoef[id][k][0];
        double x1 = (x0>0.0) ? crx.legcoef[id][k][j]/x0/(2.0*j+1.0) : 0.0;
        cout << setprecision(3) << setw(11) << lowfilter(x1);
      }
      nl();
    }
    nl();
    nl();
  }
}


/**********************************************************/
/*      Particle Emission Spectra                         */
/**********************************************************/
void outSpectrum(const int ctl, double **spc, Nucleus *n)
{
  int kmax = (int)(n->max_energy/n->de) +1;

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if (!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    int k = (int)(ncl[p].max_energy/ncl[p].de);
    if(k > kmax) kmax = k;
  }

  double gw = parmGetValue(parmBROD);
  if(gw > 0.0) kmax += (int)(3.0*gw/ncl[0].de);

  if     (ctl == 0) outSectionHead("SPECTRA FROM COMPOUND NUCLEUS");
  else if(ctl == 1) outSectionHead("PRECOMPOUND PARTICLE SPECTRA");
  else if(ctl == 2) outSectionHead("TOTAL PARTICLE EMISSION SPECTRA");

  cout << cline;
  if(ctl == 0) outZA(&n->za);
  cout << endl;

  cout << "# Emin[MeV]  Emax[MeV]";
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    cout << setw(11) << particle_name[j];
  }
  cout << endl;

  for(int k=0 ; k<=kmax ; k++){
    if(k == 0){
      outVal(11,4,0.0);
      outVal(11,4,0.5*n->de);
    }
    else{
      outVal(11,4,((double)k-0.5)*n->de);
      outVal(11,4,((double)k+0.5)*n->de);
    }
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(!n->cdt[j].status) continue;
      outVal(lowfilter(spc[j][k]));
    }
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Sum All Spectrum                                  */
/**********************************************************/
void outSpectrumSum(const int ctl, double ** spc, Nucleus * n)
{
  double sum[MAX_CHANNEL],ave[MAX_CHANNEL];
  int    kmax = (int)(n->max_energy/n->de) +1;

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    sum[j] = ave[j] = 0.0;
    if (!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    int k = (int)(ncl[p].max_energy/ncl[p].de);
    if(k > kmax) kmax = k;
  }

  for(int k=0 ; k<kmax ; k++){
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(n->cdt[j].status){
        sum[j] += spc[j][k] * n->de;
        ave[j] += spc[j][k] * n->de * k*n->de;
      }
    }
  }

  if(ctl == 0)       outSectionHead("SUM SPECTRA FOR COMPOUND");
  else if(ctl == 1)  outSectionHead("SUM OF PRECOMPOUND EMISSION");
  else if(ctl == 2)  outSectionHead("TOTAL PARTICLE EMISSION SPECTRA");

  cout << cline;
  if(ctl == 2) cout << blank;
  else         outZA(&n->za);

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    cout << setw(11) << particle_name[j];
  }
  nl();

  cout << " SpcSum[mb]" << blank;
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    outVal(lowfilter(sum[j]));
  }
  nl();

  cout << " Emean[MeV]" << blank;
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    if(sum[j] > 0.0) ave[j] = ave[j]/sum[j];
    else             ave[j] = 0.0;
    outVal(lowfilter(ave[j]));
  }
  nl();

  if(ctl == 1){
    cout << cline << endl;
    cout << cline << " TotalPreEq   Compound   Fraction" << endl;
    cout << "  Sigma[mb]";
    outVal(crx.preeq); outVal(crx.reaction); outVal(crx.preeq/crx.reaction); nl();
  }
}


/**********************************************************/
/*      Primary Gamma-Ray Spectrum                        */
/**********************************************************/
void outPrimaryGammaSpectrum(const int ng, double **spc, Nucleus *n)
{
  outSectionHead("PRIMARY GAMMA-RAY SPECTRUM FOR CAPTURE");

  cout << cline; outZA(&n->za);  cout << endl;
  cout << "#            Egam[MeV]  Sigma[mb]" << endl;

  for(int i=ng-1 ; i>=0 ; i--){
    cout << blank;
    outVal(          spc[0][i] );
    outVal(lowfilter(spc[1][i]));
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Discrete Level Population for Binary Reaction     */
/**********************************************************/
void outDiscreteLevelPopulation(double extot, Nucleus *n)
{
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    if (!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    if(ncl[p].ndisc == 0) continue;

    outSectionHead("DISCRETE LEVEL POPULATION");
    cout << cline; outZA(&ncl[p].za); cout << endl;
    cout << "#   Ex[MeV]     JP      Ecms[MeV]  Prod.[mb]" << endl;

    double sum = 0.0;
    for(int i0=0 ; i0<ncl[p].ndisc ; i0++){
      double e = extot - ncl[p].lev[i0].energy - n->cdt[j].binding_energy;
      if(e < 0.0) continue;

      outVal(11,4,ncl[p].lev[i0].energy);
      outVal(7, 1,ncl[p].lev[i0].spin);
      char p0 = (ncl[p].lev[i0].parity < 0) ? '-' : '+';
      cout << p0 << "   ";
      outVal(e);
      outVal(lowfilter(ncl[p].lpop[i0])); nl();

      sum += ncl[p].lpop[i0];
    }
    cout << blank << blank << "        sum";
    outVal(lowfilter(sum)); nl();
  }
}


/**********************************************************/
/*      Population of Continuum in Compound Nucleus       */
/**********************************************************/
void outPopulation(Nucleus *n)
{
  const int jmax = 6;

  if(n->ncont == 0) return;

  outSectionHead("CONTINUUM POPULATION");
  cout << cline; outZA(&n->za);  cout << endl;
  cout << "#   Ex[MeV]";
  for(int j=0 ; j<jmax ; j++) cout << "     "<< setw(3) << j << "   ";
  cout << "   Total   " << endl;

  for(int k=0 ; k<n->ncont ; k++){
    outVal(n->excitation[k]);

    double sum = 0.0;
    for(int j=0 ; j<=n->jmax ; j++) sum += n->pop[k][j].even+n->pop[k][j].odd;
    for(int j=0 ; j<jmax ; j++){
      double x = n->pop[k][j].even+n->pop[k][j].odd;
      outVal(lowfilter(x));
    }
    outVal(lowfilter(sum)); nl();
  }
}


/**********************************************************/
/*      Gamma Cascade                                     */
/**********************************************************/
void outGammaCascade(const int c0, Nucleus *n)
{
  outSectionHead("GAMMAS FROM DISCRETE TRANSITION");
  cout << cline;
  for(int j=1 ; j<MAX_CHANNEL ; j++) cout << setw(2) << (int)crx.prod[c0].par[j];
  cout << blank;
  outZA(&n->za); nl();
  cout << "#   Ex[MeV]     Jpi    FinalState";
  cout << "  Branching Production    Egamma Production   T(1/2)[s]" << endl;

  double tisom = parmGetValue(parmISOM);
  bool   gcut  = (tisom > 0.0) ? true : false;

  for(int i0=n->ndisc-1 ; i0>=0 ; i0--){
    if(n->lpop[i0] < output_eps && i0 > 0) continue;

    cout << setw(3) << i0;
    outVal(8,4,n->lev[i0].energy);
    outVal(7,1,n->lev[i0].spin);
    char p0 = (n->lev[i0].parity < 0) ? '-' : '+';
    cout << p0 << "   " << blank << blank;
    outVal(lowfilter(n->lpop[i0]));
    if((n->lev[i0].halflife > 0.0) && (lowfilter(n->lpop[i0]) > 0.0)){
      cout << blank <<  blank;
      outVal(n->lev[i0].halflife);
    }
    nl();

    if(i0 == 0) continue;
    if( gcut && (n->lev[i0].halflife > tisom) ) continue;
    for(int j=0 ; j<n->lev[i0].ngamma ; j++){
      int    i1 = n->lev[i0].fstate[j];
      double eg = n->lev[i0].energy-n->lev[i1].energy;
      double sg = n->lev[i0].branch[j]*n->lpop[i0];
      if(opt.internalconversion) sg *= n->lev[i0].gratio[j];

      cout << blank << blank << setw(3) << i1;
      outVal(8,4,n->lev[i1].energy);
      outVal(11,4,n->lev[i0].branch[j]);
      cout << "          ";
      outVal(11,4,eg);
      outVal(sg); nl();
    }
  }
}


/**********************************************************/
/*      Isomer Production                                 */
/**********************************************************/
void outIsomerProduction(Nucleus *n)
{
  bool   found = false;
  for(int i0=1 ; i0<n->ndisc ; i0++){
    if(n->lev[i0].halflife > thalfmin && lowfilter(n->lpop[i0]) > 0.0) found = true;
  }
  if(!found) return;

  outSectionHead("ISOMERIC STATE PRODUCTION");
  cout << cline; outZA(&n->za); cout << endl;
  cout << "#   Ex[MeV]     Jpi     T(1/2)[s]  Prod.[mb]" << endl;

  for(int i0=1 ; i0<n->ndisc ; i0++){
    if(n->lev[i0].halflife <= thalfmin) continue;

    cout << setw(3) << i0;
    outVal(8,4,n->lev[i0].energy);
    outVal(7,1,n->lev[i0].spin);
    char p0 = (n->lev[i0].parity < 0) ? '-' : '+';
    cout << p0 << "   ";
    outVal(n->lev[i0].halflife);
    outVal(lowfilter(n->lpop[i0])); nl();
  }
}


/**********************************************************/
/*      DSD Capture Cross Section                         */
/**********************************************************/
void outDSD(double sigt, double sigd, double sigc, Dcapt *dsd)
{
  outSectionHead("DSD CROSS SECTION");
  cout << cline << "       Real       Imag   NonLocal" << endl;
  cout << "    V1[MeV]";
  outVal(11,3, dsd->v1); outVal(11,3, dsd->w1); outVal(11,5, dsd->nonloc); nl();
  cout << cline << endl;
  cout << cline << "      Total     Direct Semidirect" << endl;
  cout << "  Sigma[mb]";
  outVal(sigt); outVal(sigd); outVal(sigc); nl();
}


/**********************************************************/
/*      Exciton Model Parameters                          */
/**********************************************************/
void outExcitonParameter(double **spd, double **m2, double **m2r)
{
  const int npmax = 5;

  outSectionHead("PRE-COMPOUND PARAMETERS");
  cout << cline << "   particle   gz[/MeV]   gn[/MeV] delta[MeV]  Vdep[MeV]" << endl;

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(ncl[0].cdt[j].status){
      cout << " PreEq     " << particle_name[j];
      outVal(spd[j][0]); outVal(spd[j][1]); outVal(spd[j][2]); outVal(spd[j][3]); nl();
    }
  }

  cout << cline << blank;
  for(int i=0 ; i<npmax ; i++) cout << setw(11) << i;
  cout << endl;
  cout <<  " M2 [/hbar]"; outVal(0); outVal(m2[0][0]); nl();

  cout <<  " M2R[/hbar]" << blank;
  outVal(m2r[0][0]); nl();

  for(int i=1 ; i<npmax ; i++){
    cout << blank << setw(11) << i;
    for(int j=0 ; j<=i ; j++) outVal(m2[i][j]);
    nl();

    cout << blank << blank;
    for(int j=0 ; j<=i ; j++) outVal(m2r[i][j]);
    nl();
  }
}


/**********************************************************/
/*      Fission Neutron Spectrum Model Parameters         */
/**********************************************************/
void outChanceFispecParameter(const int m, FChance *fc, Nucleus *n)
{
  outSectionHead("MULTICHANCE FISSION NEUTRON MODEL PARAMETERS");

  cout << "#";
  switch(m){
  case  0: cout << " 1st"; break;
  case  1: cout << " 2nd"; break;
  case  2: cout << " 3rd"; break;
  default: cout << " "<<setw(1) << m << "th"; break;
  }

  cout <<"          LightFF    HeavyFF   Compound" << endl;
  cout << blank;
  outZA(&fc->lf.za);
  outZA(&fc->hf.za);
  outZA(&n->za);
  nl();

  cout << "  LDP[/MeV]";
  outVal(11,4,fc->lf.a);
  outVal(11,4,fc->hf.a);
  outVal(11,4,fc->ac);
  nl();

  cout << "  Tmax[MeV]";
  outVal(11,4,fc->lf.tmax);
  outVal(11,4,fc->hf.tmax);
  outVal(11,4,fc->tmax);
  nl();

  cout << cline << endl;
  cout << cline << " Efiss[MeV]   TKE[MeV]   TXE[MeV]  Egam[MeV] Exfis[MeV]" << endl;
  cout << "     Energy";
  outVal(fc->etotal); outVal(fc->tke); outVal(fc->txe); outVal(fc->egamma); outVal(fc->exfiss);
  nl();

  cout << cline << endl;
  cout << cline << "         Rt   Nu-Ratio     tpdf_s anisotropy" << endl;
  cout << "  Parameter";
  outVal(11,4,fc->rt); outVal(11,4,fc->nuratio);
  outVal(11,4,fc->tps); outVal(11,4,fc->anisotropy);
  nl();
}


/**********************************************************/
/*      Multi-Chance Fission Neutron Spectrum             */
/**********************************************************/
void outChanceFispec(const int m, double *e, double **s, FChance *fc, Nucleus *n)
{
  outSectionHead("MULTICHANCE FISSION NEUTRON SPECTRUM");
  cout << cline; outZA(&fc->lf.za); outZA(&fc->hf.za);
  cout << blank; outZA(&n->za); nl();
  cout << "#    E[MeV]    LightFF    HeavyFF    Average   Compound      Total" << endl;

  for(int k=0 ; k<m ; k++){
    outVal(e[k]);
    for(int i=0 ; i<5 ; i++) outVal(lowfilter(s[i][k]));
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Multi-Chance Fission Neutron Average Energy       */
/**********************************************************/
void outChanceFispecEnergy(double *ecms)
{
  outSectionHead("MULTICHANCE FISSION NEUTRON ENERGIES");

  cout << cline << "    LightFF    HeavyFF    Average   Compound      Total" << endl;
  cout << "  Ecms[MeV]";
  for(int i=0 ; i<5 ; i++) outVal(11,4,ecms[i]);
  nl();
}


/**********************************************************/
/*      Total Fission Neutron Energy (<E> and Nu-bar)     */
/**********************************************************/
void outFissionNeutronEnergy(const int mc, double e, double nu, double *pf, FNSpec *fns)
{
  outSectionHead("AVERAGE FISSION NEUTRON ENERGY AND MULTIPLICITY");

  cout << cline << "  FissProb. Epref[MeV]     Nu-bar Efiss[MeV]   TKE[MeV]   TXE[MeV]  Egam[MeV] Exfis[MeV]" << endl;

  const int neave = 5;
  double eave[neave];
  for(int j=0 ; j<neave ; j++) eave[j] = 0.0;

  for(int i=0 ; i<=mc ; i++){
    cout << " FisChance" << setw(1) <<i+1;
    outVal(pf[i]);
    outVal(11,4,fns->fc[i].ecms[3]);
    outVal(11,4,fns->fc[i].nubar + i);
    outVal(11,4,fns->fc[i].etotal);
    outVal(11,4,fns->fc[i].tke);
    outVal(11,4,fns->fc[i].txe);
    outVal(11,4,fns->fc[i].egamma);
    outVal(11,4,fns->fc[i].exfiss);
    nl();

    int j = 0;
    eave[j++] += pf[i] * fns->fc[i].etotal;
    eave[j++] += pf[i] * fns->fc[i].tke;
    eave[j++] += pf[i] * fns->fc[i].txe;
    eave[j++] += pf[i] * fns->fc[i].egamma;
    eave[j]   += pf[i] * fns->fc[i].exfiss;
  }

  cout << "    Average" << blank;
  outVal(11,4,e);
  outVal(11,4,nu);
  for(int j=0 ; j<neave ; j++) outVal(11,4,eave[j]);
  nl();

/*
  cout << " FissEnergy";
  outVal(labE);
  for(int i=0 ; i<4 ; i ++){
    if(i <= mc) outVal(pf[i]); else{ outVal(11,4,0.0); }
  }
  for(int i=0 ; i<4 ; i ++){
    if(i <= mc) outVal(fns->fc[i].tke); else{ outVal(11,4,0.0); }
  }
  for(int i=0 ; i<4 ; i ++){
    if(i <= mc) outVal(fns->fc[i].exfiss); else{ outVal(11,4,0.0); }
  }
  for(int i=0 ; i<4 ; i ++){
    if(i <= mc) outVal(fns->fc[i].ecms[3]); else{ outVal(11,4,0.0); }
  }
  nl();
*/
}


/**********************************************************/
/*      Total Fission Neutron Spectrum (Chi)              */
/**********************************************************/
void outFissionNeutronSpectrum(const int n, double *e, double *s, double t)
{
  outSectionHead("FISSION NEUTRON SPECTRUM");

  if(t > 0.0){
    cout << "#    E[MeV]  Chi[/MeV]   /Maxwell" << endl;
    double c = 2.0/sqrt(PI * t*t*t);
    for(int k=0 ; k<n ; k++){
      outVal(e[k]);
      outVal(lowfilter(s[k]));
      double y = c * exp(-e[k]/t) * sqrt(e[k]);
      outVal(lowfilter(s[k] / y));
      nl();
    }
  }
  else{
    cout << "#    E[MeV]  Chi[/MeV]" << endl;
    for(int k=0 ; k<n ; k++){
      outVal(e[k]);
      outVal(lowfilter(s[k]));
      nl();
    }
  }
  nl();
  nl();
}


/**********************************************************/
/*      Adjustable Model Parameters                       */
/**********************************************************/
void outParameter()
{
  if(adj.nparm() <= 0) return;

  int k = 0;
  for(int i=0 ; i<adj.nparm() ; i++){
    if(adj.parm[i].type == parmLCUT) continue;
    if(adj.parm[i].type == parmGSTR) continue;
    if(adj.parm[i].type == parmISOM) continue;
    if(adj.parm[i].type == parmD0SR) continue;
    if(adj.parm[i].type == parmESET) continue;
    k++;
  }
  if(k == 0) return;

  outSectionHead("ADJUSTED PARAMETERS");
  cout << "# Parameter    Nucleus   Particle     Factor" << endl;

  for(int i=0 ; i<adj.nparm() ; i++){
    if(adj.parm[i].type == parmLCUT) continue;
    if(adj.parm[i].type == parmGSTR) continue;
    if(adj.parm[i].type == parmISOM) continue;
    if(adj.parm[i].type == parmD0SR) continue;
    if(adj.parm[i].type == parmESET) continue;

    int      t  = adj.parm[i].type;
    int      p  = adj.parm[i].particle; if(p > 8) p = 8;
    ZAnumber za = adj.parm[i].za;
    double   f  = adj.parm[i].getfactor();

    cout << "    " << pname[t] << "   ";
    switch(t){
    case parmM2  :
    case parmM2R :
    case parmKO  :
    case parmDSDV:
    case parmGDRM1:
                    cout << blank << blank; 
                    break;
    case parmOV  :
    case parmOW  :
    case parmORV :
    case parmORW :
    case parmOAV :
    case parmOAW :
    case parmOC  :
    case parmSD  :  cout << blank << particle_name[p];
                    break;
    case parmLD  :
    case parmSPIN:
    case parmPAIR:
    case parmPDST:
    case parmFL1 :
    case parmFL2 :
    case parmFL3 :
                    outZA(&za);   cout << blank;
                    break;
    case parmTJ  :  outZA(&za);   cout << particle_name[p];
                    break;
    default      :  break;
    }
    outVal(f); nl();
  }
}

