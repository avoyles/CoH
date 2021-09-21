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

using namespace std;

#include "beoh.h"
#include "beohoutput.h"
#include "nucleus.h"
#include "elements.h"
#include "outformat.h"


/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/

extern string version;

static string particle_name[7]={
   "      gamma","    neutron","     proton","      alpha",
   "   deuteron","     triton","     helion"};

static string p_name = "gnpadth";
static double Qbeta = 0.0;
static ZAnumber targZA(0,0);


/**********************************************************/
/*      Print Banner                                      */
/**********************************************************/
void outBanner()
{
  outSectionHead(&version[0]);
  cout
    <<"#    oooooooooo.                     ooooo   ooooo\n"
    <<"#     888     Y8b                     888     888      In spring, it's the dawn.\n"
    <<"#     888     d8P  .oooo.   .ooooo.   888     888      The mountain rim is slowly getting light,\n"
    <<"#     888ooooooB  d88  88b d88   88b  888ooooo888      and thin purple cloud wisps there.\n"
    <<"#     888     Y8b 888oo888 888   888  888     888\n"
    <<"#     888     d8P 888      888   888  888     888                                The Pillow Book\n"
    <<"#    o8888bood8P   Y8bodP'  Y8boodP' o888o   o888o                                  Sei Shonagon\n"
    <<"#\n";
}


/**********************************************************/
/*      Header                                            */
/**********************************************************/
void outTitle(char *str)
{
  outBanner();
  outSectionHead(&version[0]);
  cout <<cline;
  cout << &str[10] << endl;
}


/**********************************************************/
/*      Write  Z and A                                    */
/**********************************************************/
void outZA(ZAnumber *za)
{
  char element[3];

  if(za->getZ() >= N_ELEMENTS){
    element[0] = element[1] = ' ';
  }
  else{
    strncpy(element,element_name[za->getZ()].c_str(),2);
    if(strlen(element_name[za->getZ()].c_str())==1) element[1] = ' ';
  }
  element[2] = '\0';

  cout << "  ";
  cout << setw(3) << setfill('0') << za->getZ() << '-';
  cout << setw(3) << setfill('0') << za->getA();
  cout << setw(2) << element;
  cout << setfill(' ');
}


/**********************************************************/
/*      Write  Header of Each Section                     */
/**********************************************************/
void outSectionHead(const char *sec)
{
  int n = strlen(sec);
  const int DisplayWidth = 99;

  string bar1="#";

  if(Qbeta == 0.0){
    for(int i=0 ; i<DisplayWidth-1 ; i++) bar1 += "#";
    cout << bar1 << endl;
  }else{
    for(int i=0 ; i<DisplayWidth-27 ; i++) bar1 += ".";
    cout << bar1 << "/";  outZA(&targZA);  cout << "  //";  outVal(10,6,Qbeta); nl();
  }

  string bar2 = "#";
  for(int i=0 ; i<DisplayWidth-n-1 ; i++) bar2 += " ";
  cout << bar2 << sec << endl;
}


/**********************************************************/
/*      System Parameters                                 */
/**********************************************************/
void outSystem(CalcMode mode, System *sys)
{
  /*** save Z, A, E for the first time */
  targZA.setZA(sys->target.getZ(),sys->target.getA());
  Qbeta  = sys->ex_total;

  outSectionHead("SYSTEM PARAMETERS");
  if(mode == betadecay){
    cout << cline << "  AtomicNum    MassNum     Q-beta" << endl;
    cout << "# Precursor"; outVal(11,sys->target.getZ()); outVal(11,sys->target.getA()); outVal(11,4,sys->ex_total); nl();
    cout << "#  Compound"; outVal(11,sys->compound.getZ()); outVal(11,sys->compound.getA()); nl();
  }
  else if(mode == statdecay){
    cout << cline << "  AtomicNum    MassNum Excitation";
    if(sys->beta2 != 0.0) cout << "      beta2";
    nl();

    cout << "#  Compound"; outVal(11,sys->compound.getZ()); outVal(11,sys->compound.getA()); outVal(11,4,sys->ex_total);
    if(sys->beta2 != 0.0) outVal(11,4,sys->beta2);
    nl();
  }
  else if(mode == fissiondecay || mode == fissionspec || mode == cumulativeyield){
    cout << cline << "  AtomicNum    MassNum Excitation" << endl;
    cout << "#  Target  "; outVal(11,sys->target.getZ()); outVal(11,sys->target.getA()); nl();
    cout << "#  Compound"; outVal(11,sys->compound.getZ()); outVal(11,sys->compound.getA()); outVal(11,4,sys->ex_total); nl();
  }
}


/**********************************************************/
/*      Parameters fof Fission Fragment Decay Mode        */
/**********************************************************/
void outFissionFragment(const int af, FFragData *fdt)
{
  outSectionHead("FISSION FRAGMENT DECAY PARAMETER");

  cout << cline << "Yield Cutoff          "; outVal(11,fdt->ycutoff); nl();
  cout << cline; nl();

  cout << "# FisChance   Fraction      Rtemp SpinFactor   Zpfac(Z)   Zpfac(N)        TKE    Ex_fiss   Eprefiss" << endl;
  for(int n=0 ; n<fdt->getFissionChance() ; n++){
    cout << "# " << setw(2) << n+1 << setw(7) << af - n;
    outVal(11,fdt->mc[n].fraction);
    outVal(11,5,fdt->mc[n].rt);
    outVal(11,5,fdt->mc[n].spinfactor);
    outVal(11,5,fdt->mc[n].ZpFactor[0]);
    outVal(11,5,fdt->mc[n].ZpFactor[1]);
    outVal(11,fdt->mc[n].tke);
    outVal(11,fdt->mc[n].exfis);
    outVal(11,fdt->mc[n].eprefis);
    nl();
  }
  nl();

  cout << "# Gaussian Parameters    fraction      width     center"; nl();
  for(int n=0 ; n<fdt->getFissionChance() ; n++){
    for(int i=0 ; i<4 ; i++){
      if(i == 0)  cout << "# " << setw(9) << n+1;
      else cout << cline;
      cout << setw(11) << i+1;
      outVal(11,5,fdt->mc[n].GaussFract[i]);
      outVal(11,5,fdt->mc[n].GaussSigma[i]);
      outVal(11,5,fdt->mc[n].GaussDelta[i]); nl();
    }
  }
  nl();
}


/**********************************************************/
/*      Target State                                      */
/**********************************************************/
void outTargetState(NuclearStructure *nst, const int isostate)
{
  outSectionHead("PRECURSOR STATE");
  cout << cline << " Excitation  Spin" << endl;

  outZA(&nst->za);
  outVal(nst->lev[isostate].energy);
  outVal(5,1,nst->lev[isostate].spin);
  char p = (nst->lev[isostate].parity < 0) ? '-' : '+';
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
    if(ncl[i].ncont>0){
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
void outGDR(GDR *gdr)
{
  outSectionHead("GIANT DIPOLE RESONANCE DATA");
  cout << cline << "Energy[MeV] Width[MeV] Sigma0[mb]" << endl;
  for(int i=0 ; i<MAX_GDR ; i++){
    if(gdr[i].getEnergy()>0.0){
      cout << "         " << gdr[i].getEM() << setw(1) << gdr[i].getL();
      outVal(11,4,gdr[i].getEnergy());
      outVal(11,4,gdr[i].getWidth());
      outVal(11,6,gdr[i].getSigma()); nl();
    }
  }
}


/**********************************************************/
/*      Fission Barriers                                  */
/**********************************************************/
void outFissionBarrier(const int n)
{
  outSectionHead("FISSION BARRIER PARAMETERS");
  cout << cline << blank << "Height[MeV] Width[MeV] KBand[MeV]" << endl;
  for(int i=0 ; i<n ; i++){
    for(int m=0 ; m<MAX_HUMP ; m++){
      if(ncl[i].fissile){
        if(ncl[i].fission->barrier[m].height > 0.0){
          outZA(&ncl[i].za);
          if     (m==0) cout << " 1stBarrier";
          else if(m==1) cout << " 2ndBarrier";
          else if(m==2) cout << " 3rdBarrier";
          else          cout << blank;
          outVal(11,4,ncl[i].fission->barrier[m].height);
          outVal(11,4,ncl[i].fission->barrier[m].curvature); nl();

          for(int k=0 ; k<ncl[i].fission->barrier[m].nband ; k++){
            cout << blank << blank << blank << blank;
            outVal(ncl[i].fission->barrier[m].kband[k].excitation);
            outVal(5,1,ncl[i].fission->barrier[m].kband[k].k2/2.0);
            char p = (ncl[i].fission->barrier[m].kband[k].parity < 0) ? '-' : '+';
            cout << p << endl;
          }
        }
      }
    }
  }
}


/**********************************************************/
/*      Beta Decay Strength                               */
/**********************************************************/
void outBeta(ZAnumber *res, Beta *gts, Beta *ens)
{
  outSectionHead("BETA DECAY FINAL STATE");
  cout << "# Daughter     GTstate ENSDFlevel" << endl;
  outZA(res);
  outVal(gts->nstate);
  outVal(ens->nstate);
  cout << endl;
}


/**********************************************************/
/*      Beta Decay Strength Distribution                  */
/**********************************************************/
void outBetaProfile(BetaProfile *bpf)
{
  outSectionHead("BETA STRENGTH DISTRIBUTION");

  cout << "# DataFrom " << bpf->source << endl;
  cout << "# Fraction    Discrete  Continuum" << endl;
  cout << blank;
  outVal(11,4,bpf->tdisc); outVal(11,4,bpf->tcont);
  cout << endl;

/*
  if(bpf->ndisc > 0){
    cout << "# Discrete" << endl;
    cout << "#     LeveL     Spin   Population" << endl;
    for(int i=0 ; i<bpf->ndisc ; i++){
      outVal(11,4,bpf->lev[i].energy);
      outVal(8,1,bpf->lev[i].spin);
      cout << ((bpf->lev[i].parity < 0) ? " - " : " + " );
      outVal(bpf->rdisc[i]);
      cout << endl;
    }
  }
  if(bpf->ncont > 0){
    cout << "# Continuum" << endl;
    cout << "#      Emin       Emax Population" << endl;
    for(int i=bpf->ncont-1 ; i>=0 ; i--){
      outVal(11,3,bpf->excitation[i+1]);
      outVal(11,3,bpf->excitation[i]);
      outVal(bpf->rcont[i]);
      cout << endl;
    }
  }
*/
}


/**********************************************************/
/*      Ground State Production Rate                      */
/**********************************************************/
void outGSProduction(const int n)
{
  outSectionHead("GROUND STATE PRODUCTION RATE");
  cout << cline << "   Residual  Prod.Rate   Particles" << endl;
  cout << cline << blank << blank << "  ";
  for(int j=1 ; j<(int)p_name.length() ; j++){
    cout << setw(2)<< p_name.substr(j,1);
  }
  nl();

  double sum=0.0;
  for(int i=0 ; i<n ; i++){
    cout << " Production";
    outZA(&ncl[i].za);
    outVal(lowfilter(crx.prod[i].xsec));
    cout << "  ";
    sum += crx.prod[i].xsec;

    for(int j=1 ; j<MAX_CHANNEL ; j++) cout << setw(2) << (int)crx.prod[i].par[j];
    nl();
  }
  cout << cline << "        sum";
  outVal(lowfilter(sum)); nl();
}


/**********************************************************/
/*      Fission Rate                                      */
/**********************************************************/
void outFission(int n)
{
  outSectionHead("FISSON RATE");
  cout << cline << "   Nucleus FissionRate Fiss.Prob." << endl;

  double sum = 0.0;
  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].fiss>0) sum += crx.prod[i].fiss;
  }

  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].fiss>0){
      cout << " FisChance" << setw(1) <<i+1;
      outZA(&ncl[i].za);
      outVal(lowfilter(crx.prod[i].fiss));
      outVal(11,3,crx.prod[i].fiss/sum*100.0);
      nl();
    }
  }
  cout << blank << " TotFission";
  outVal(lowfilter(sum)); nl();
}


/**********************************************************/
/*      Particle Emission Spectra                         */
/**********************************************************/
void outSpectrum(const bool betacalc, const double de, double **spc)
{
  int k0 = beohZeroCut(spc);

  outSectionHead("TOTAL PARTICLE EMISSION SPECTRA");
  cout << "#      Emin       Emax"
       << "  Spectra / Decay / MeV" << endl;
  cout << "#     [MeV]      [MeV]";

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    cout << setw(11) << particle_name[j];
  }
  if(betacalc) cout << "   electron   neutrino";
  nl();

  for(int k=0 ; k<=k0 ; k++){
    if(k==0){
      outVal(11,4,0.0);
      outVal(11,4,0.5*de);
    }
    else{
      outVal(11,4,((double)k-0.5)*de);
      outVal(11,4,((double)k+0.5)*de);
    }

    for(int j=0 ; j<MAX_CHANNEL ; j++){
      outVal(lowfilter(spc[j][k]));
    }
    if(betacalc){
      outVal(lowfilter(spc[MAX_CHANNEL  ][k]));
      outVal(lowfilter(spc[MAX_CHANNEL+1][k]));
    }
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Sum of Spectra and Average Energies               */
/**********************************************************/
void outSpectrumSum(const bool betacalc, const double de, double **spc)
{
  double et[MAX_CHANNEL+2], ea[MAX_CHANNEL+2], em[MAX_CHANNEL+2];

  int k0 = beohZeroCut(spc);
  double e, e0, e1;

  int cm = MAX_CHANNEL;
  if(betacalc) cm += 2;

  for(int c=0 ; c<cm ; c++){
    em[c] = et[c] = 0.0;
    for(int k=0 ; k<=k0 ; k++){

      e0 = (k > 0)  ? (k - 0.5) * de : 0;
      e1 = (k + 0.5) * de;
      e  = (e0 + e1) * 0.5;

      em[c] += spc[c][k] * de;
      et[c] += spc[c][k] * de * e;
    }
    ea[c] = (em[c]>0.0) ? et[c]/em[c] : 0.0;
  }


  outSectionHead("SUM SPECTRA AND AVERAGE ENERGIES");

  cout << cline << blank;
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    cout << setw(11) << particle_name[j];
  }
  if(betacalc) cout << "   Electron   Neutrino";
  nl();

  cout << "      TotalEnergy[MeV]";
  double etot = 0.0;
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    outVal(lowfilter(et[j]));
    etot += et[j];
  }
  if(betacalc){
    outVal(lowfilter(et[MAX_CHANNEL  ]));  etot += et[MAX_CHANNEL  ];
    outVal(lowfilter(et[MAX_CHANNEL+1]));  etot += et[MAX_CHANNEL+1];
  }
  nl();

  cout << "    AverageEnergy[MeV]";
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    outVal(lowfilter(ea[j]));
  }
  if(betacalc){
    outVal(lowfilter(ea[MAX_CHANNEL  ]));
    outVal(lowfilter(ea[MAX_CHANNEL+1]));
  }
  nl();

  cout << "          Multiplicity";
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    outVal(lowfilter(em[j]));
  }
  if(betacalc){
    outVal(lowfilter(em[MAX_CHANNEL  ]));
    outVal(lowfilter(em[MAX_CHANNEL+1]));
  }
  nl();
  nl();

  cout << "    EnergyRelease[MeV]";
  outVal(etot); nl();
  nl();
  nl();
}


/**********************************************************/
/*      Neutron Emission Spectra in Lab Frame             */
/**********************************************************/
void outSpectrumLab(const int k0, const double de, double *spl)
{
  outSectionHead("NEUTRON EMISSION SPECTRA IN LAB FRAME");
  cout << "#      Emin       Emax"
       << "  Spectra / Decay / MeV" << endl;
  cout << "#     [MeV]      [MeV]";

  cout << setw(11) << particle_name[1];
  nl();

  for(int k=0 ; k<=k0 ; k++){
    if(k==0){
      outVal(11,4,0.0);
      outVal(11,4,0.5*de);
    }
    else{
      outVal(11,4,((double)k-0.5)*de);
      outVal(11,4,((double)k+0.5)*de);
    }

    outVal(lowfilter(spl[k]));
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Gamma Cascade                                     */
/**********************************************************/
void outGammaCascade(Nucleus *n)
{
  outSectionHead("GAMMAS FROM DISCRETE TRANSITION");
  cout << cline; outZA(&n->za); cout << endl;
  cout << "#   Ex[MeV]     Jpi    FinalState";
  cout << "  Branching Production    Egamma Production   T(1/2)[s]" << endl;

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
    for(int j=0 ; j<n->lev[i0].ngamma ; j++){
      int    i1 = n->lev[i0].fstate[j];
      double eg = n->lev[i0].energy-n->lev[i1].energy;
      double sg = n->lev[i0].branch[j]*n->lpop[i0];

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
/*      Total Residual Nucleus Production Cross Section   */
/**********************************************************/
void outTotalResidual(CumulativeResidualProduct *res)
{
  outSectionHead("STABLE AND LONG-LIVED STATE PRODUCTIONS");
  if(res->getNcurrent() == 0) return;

  cout << "#       Residual    Ex[MeV]     Jpi     T(1/2)[s] Prod.Ratio Meta" << endl;

  int c = 1;
  double sum0 = 0.0;
  double sum1 = 0.0;
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

        if(res->rp[i].metaflag == 0) sum0 += res->rp[i].production;
        sum1 += res->rp[i].production;
      }
    }
  }
  cout << cline << blank << blank << "         g.s.sum";
  outVal(lowfilter(sum0)); nl();
  cout << cline << blank << blank << "             sum";
  outVal(lowfilter(sum1)); nl();
  nl(); nl();
}


/**********************************************************/
/*      Isomeric Ratios for Long-Lived Nuclides           */
/**********************************************************/
void outIsomericRatio(CumulativeResidualProduct *res)
{
  outSectionHead("ISOMERIC RATIOS");
  if(res->getNcurrent() == 0) return;

  cout << "#       Residual    Ex[MeV]     Jpi     T(1/2)[s] m/(g+m)    m/g" << endl;

  int c = 1;
  for(unsigned int z=res->zmin ; z<=res->zmax ; z++){
    for(unsigned int a=res->amin ; a<=res->amax ; a++){
      ZAnumber za(z,a);

      /*** look for ground state of the same ZA */
      for(int i0=0 ; i0<res->getNcurrent() ; i0++){
        if(res->rp[i0].za != za) continue;
        if(res->rp[i0].metaflag != 0) continue;

        double pg = res->rp[i0].production;

        /*** sum g.s. and isomers */
        double pt = res->rp[i0].production;
        for(int i1=0 ; i1<res->getNcurrent() ; i1++){
          if((res->rp[i1].za == za) && (res->rp[i1].metaflag > 0)){
            pt += res->rp[i1].production;
          }
        }

        /*** look at the isomer with the same ZA */
        for(int i1=0 ; i1<res->getNcurrent() ; i1++){
          if((res->rp[i1].za != za) || (res->rp[i1].metaflag == 0)) continue;

          double r0 = (pt == 0.0) ? 0.0 : res->rp[i1].production / pt;
          double r1 = (pg == 0.0) ? 0.0 : res->rp[i1].production / pg;

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
  nl(); nl();
}


/**********************************************************/
/*      Print Fission Product Yield Data                  */
/**********************************************************/
void outFissionProductYield(const int nuktotal, Isotope *nuk)
{
  outSectionHead("INDEPENDENT AND CUMULATIVE FISSION YIELDS");
  if(nuktotal == 0) return;

  /*** find min/max Z and A */
  int z0 = 100, z1 = 0, a0 = 1000, a1 = 0;
  for(int k=0 ; k<nuktotal ; k++){
    int z = nuk[k].getZ();
    int a = nuk[k].getA();
    if(z < z0) z0 = z;
    if(z > z1) z1 = z;
    if(a < a0) a0 = a;
    if(a > a1) a1 = a;
  }

  cout << "#        Isotope M     Indep.FPY     Cumul.FPY  T(1/2)[s]" << endl;

  int c = 1;
  double s0 = 0.0, s1 = 0.0;
  for(int z = z0 ; z <= z1 ; z++){
    for(int a = a0 ; a <= a1 ; a++){

      ZAnumber za(z,a);

      for(int k=0 ; k<nuktotal ; k++){
        if((unsigned int)z != nuk[k].getZ() || (unsigned int)a != nuk[k].getA()) continue;
        if((nuk[k].initial == 0.0) && (nuk[k].yield == 0.0)) continue;

        double t2 = (nuk[k].getLambda() == 0.0) ? 0.0 : log(2.0)/nuk[k].getLambda();

        outVal(5,c++);
        outZA(&za);
        cout << setw(2) << nuk[k].getM();
        outVal(14,nuk[k].initial);
        outVal(14,nuk[k].yield);
        outVal(11,t2);
        nl();
        s0 += nuk[k].initial;
        s1 += nuk[k].yield;
      }
    }
  }
  nl();

  cout << "#            Total";
  outVal(14,s0);
  outVal(14,s1);
  nl(); nl();
}


/**********************************************************/
/*      Calculate and Print Mass and Chain Yield          */
/**********************************************************/
void outFissionProductChainYield(const int nuktotal, Isotope *nuk)
{
  outSectionHead("MASS AND CHAIN FISSION YIELDS");
  if(nuktotal == 0) return;

  /*** find min/max A */
  int a0 = 1000, a1 = 0;
  for(int k=0 ; k<nuktotal ; k++){
    int a = nuk[k].getA();
    if(a < a0) a0 = a;
    if(a > a1) a1 = a;
  }

  cout <<"#        A    IndepYield     MassYield    ChainYield" << endl;

  double s0, s1, s2, c0, c1, c2;
  /*** nuclide which has the longest half-life in a mass chain */
  s0 = s1 = s2 = 0.0;
  for(int a=a0 ; a<=a1 ; a++){

    c0 = c1 = c2 = 0.0;
    for(int k=0 ; k<nuktotal ; k++){
      if(nuk[k].getA() != (unsigned int) a) continue;

      /*** initial yields */
      c0 += nuk[k].initial;

      /*** just add all the yields of the same mass */
      c1 += nuk[k].yield;

      /*** find a stable nuclide in a decay chain */
      if(nuk[k].isstable()) c2 += nuk[k].yield;
    }

    cout << setw(10) << a;
    outVal(14,c0);
    outVal(14,c1);
    outVal(14,c2); nl();

    s0 += c0;
    s1 += c1;
    s2 += c2;
  }
  nl();

  cout << "#    Total";
  outVal(14,s0);
  outVal(14,s1);
  outVal(14,s2);
  nl(); nl();
}


/**********************************************************/
/*      Delayed Neutron Yield                             */
/**********************************************************/
void outDelayedNeutronYield(const int nuktotal, Isotope *nuk)
{
  outSectionHead("DELAYED NEUTRON YIELD");
  if(nuktotal == 0) return;

  /*** find min/max Z and A */
  int z0 = 100, z1 = 0, a0 = 1000, a1 = 0;
  for(int k=0 ; k<nuktotal ; k++){
    int z = nuk[k].getZ();
    int a = nuk[k].getA();
    if(z < z0) z0 = z;
    if(z > z1) z1 = z;
    if(a < a0) a0 = a;
    if(a > a1) a1 = a;
  }

  cout << "#     DN-Emitter M  T(1/2)[s]  DN-Branch     Cumul.FPY      DN-Yield" << endl;

  int c = 1;
  double dn = 0.0, dgroup6[6], lambda6[6];
  for(int i=0 ; i<6 ; i++) dgroup6[i] = lambda6[i] = 0.0;

  for(int z = z0 ; z <= z1 ; z++){
    for(int a = a0 ; a <= a1 ; a++){

      ZAnumber za(z,a);

      for(int k=0 ; k<nuktotal ; k++){
        if(nuk[k].isstable()) continue;
        if((unsigned int)z != nuk[k].getZ() || (unsigned int)a != nuk[k].getA()) continue;

        double t2 = (nuk[k].getLambda() == 0.0) ? 0.0 : log(2.0)/nuk[k].getLambda();

        /*** search for delayed neutron emitters */
        for(int m=0 ; m<nuk[k].mode.getNDM() ; m++){
          int j = nuk[k].mode.getNext(m);
          /*** beta-decay and delayed neutron emission, Zi+1 = Zj and Ai != Aj  */
          if( (nuk[k].getZ()+1 == nuk[j].getZ()) && (nuk[k].getA() != nuk[j].getA()) ){

            /*** (cumulative yield) x (number of neutrons) x (branching ratio) */
            double x = nuk[k].yield * (nuk[k].getA() - nuk[j].getA()) * nuk[k].mode.getBranch(m);

            outVal(5,c++);
            outZA(&za);
            cout << setw(2) << nuk[k].getM();
            outVal(11,t2);
            outVal(11,nuk[k].mode.getBranch(m));
            outVal(14,nuk[k].yield);
            outVal(14,x);
            nl();

            dn += x;

            /*** six group constant */
            int gid = 0;
            if(     t2 > 40.0) gid = 0;
            else if(t2 >  8.0) gid = 1;
            else if(t2 >  3.0) gid = 2;
            else if(t2 >  1.0) gid = 3;
            else if(t2 >  0.3) gid = 4;
            else               gid = 5;

            dgroup6[gid] += x;
            lambda6[gid] += x * nuk[k].getLambda();
          }
        }
      }
    }
  }

  /*** decay constant weighted average */
  for(int i=0 ; i<6 ; i++){
    if(dgroup6[i] > 0.0) lambda6[i] = lambda6[i] / dgroup6[i];
    else lambda6[i] = 0.0;
  }

  nl();
  cout << "#                          Total Delayed Neutron Yield";
  outVal(14,dn);
  nl();
  nl();

  cout << "#         DN-Group DecayConstant   Group-Yield" << endl;
  for(int i=0 ; i<6 ; i++){
    cout << setw(18) << i+1;
    outVal(14,lambda6[i]);
    outVal(14,dgroup6[i]);
    nl();
  }
  nl(); nl();
}

