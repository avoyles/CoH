/******************************************************************************/
/*  HObasis.cpp                                                               */
/*        Harmonic Oscillator basis for cylindrical coordinate                */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "HObasis.h"

static void   HOSphericalBasis (SphericalBasis *);
static double OverlapFunction  (const int, const int, const int,
                                Basis *, SphericalBasis *, Operator *);
static double OverlapFunctionA (int, int, double);
static double OverlapFunctionB (int, int, int, int, int);

extern double *fact;

/**********************************************************/
/*      Cylindrical Harmonic Oscillator Basis             */
/*      --------                                          */
/*      generate single particle states                   */
/*      based on the cylindrical HO                       */
/**********************************************************/
int HOCylindricalBasis(Basis *bs)
{
  /*** scaling constants */
  bs->bp  = bs->b * pow(bs->q, 1.0/6.0);
  bs->bz  = bs->b / pow(bs->q, 1.0/3.0);
  bs->bp2 = bs->bp * bs->bp;
  bs->bz2 = bs->bz * bs->bz;

  double q13 = pow(bs->q, 1.0/3.0);
  double q23 = q13*q13;
  double x   = (bs->getN0() + 2.0)/q13 - 1.0 - 0.5/bs->q;

  bs->nzmax = (int)(pow(bs->q, 2.0/3.0) * (bs->getN0()+2.0) - bs->q - 0.5);
  bs->nrmax = (int)(0.5 * x);
  bs->lambdamax = (int)x;

  double e0 = bs->getN0() + 2.0;
  int omega_max = 4*bs->nrmax+3;

  /*** count number of total vectors and blocks (same omega) */
  bs->nv = 0;
  bs->nb = 0;
  for(int m=1 ; m<=omega_max ; m+=2){
    int l = (m+1)/2;
    bool found = false;
    for(int is=-1 ; is<=1 ; is+=2){
      for(int ir=0 ; ir<=bs->nrmax ; ir++){
        for(int iz=0 ; iz<=bs->nzmax ; iz++){
          for(int il=l-1 ; il<=l ; il++){
            if((2*il+is) == m){
              double e1 = (2.0*ir + il + 1.0)*q13 + (iz + 0.5)/q23;
              if(e1 <= e0){
                bs->nv ++;
                found = true;
              }
            }
          }
        }
      }
    }
    if(found) bs->nb ++;
  }

  if( (bs->nv > 0) && (bs->nb > 0) ) bs->memalloc();
  else return(-1);

  int i  = 0;
  int nc = 0;
  bs->index[0] = 0;
  for(int m=1 ; m<=omega_max ; m+=2){
    int l = (m+1)/2;
    int k = 0;
    for(int is=-1 ; is<=1 ; is+=2){
      for(int ir=0 ; ir<=bs->nrmax ; ir++){
        for(int iz=0 ; iz<=bs->nzmax ; iz++){
          for(int il=l-1 ; il<=l ; il++){
            if((2*il+is) == m){
              double e1 = (2.0*ir + il + 1.0)*q13 + (iz + 0.5)/q23;
              if(e1 <= e0){
                int p = ((iz + 2*ir + il)%2 == 0) ? 1 : -1;
                bs->state[i].set(ir,iz,il,is,m,p);
                i++;
                k++;
              }
            }
          }
        }
      }
    }
    if(k > 0){
      bs->omega[nc] = m;                   // omega of the block
      bs->ndata[nc] = k;                   // number of states in the block
      bs->index[nc+1] = bs->index[nc] + k; // block starting number
      nc++;
    }
  }

/*
  cout << "# Basis Parameters"<< endl;
  cout << "#   N0 = " << setw(2) << bs->getN0() << endl;
  cout << "#   b  = " << setw(12) << bs->b << endl;
  cout << "#   q  = " << setw(12) << bs->q << endl;

  cout << "#  Nr max = " << setw(12) << bs->nrmax << endl;
  cout << "#  Nz max = " << setw(12) << bs->nzmax << endl;
  cout << "#  L  max = " << setw(12) << bs->lambdamax << endl;

  cout << "# Basis States       " << setw(4) << bs->nv << endl;
  cout << "# Diagonal Blocks    " << setw(4) << bs->nb << endl;
  cout << "# Size of Each Block ";
  for(int k=0 ; k<bs->nb ; k++) cout << setw(4) << bs->ndata[k];
  cout << endl;
  cout << "# Max Jz Value       " << setw(4) << bs->state[bs->nv-1].getJ() << "/2" << endl;


  i = 0;
  for(int k=0 ; k<bs->nb ; k++){
    cout << setw(4) << k;
    cout << "  N " << setw(4) << bs->ndata[k];
    cout << "  Omega " << setw(3) << bs->omega[k]<<"/2";
    cout << "  Index " << setw(4) << bs->index[k] << setw(4) << bs->index[k+1];
    cout << endl;

    for(int j=0 ; j<bs->ndata[k] ; j++){
      cout << " " << setw(4) << i;
      cout << " : nz" << setw(3) << bs->state[i].getNz();
      cout << " : nr" << setw(3) << bs->state[i].getNr();
      cout << " : ";
      cout << setw(3) << bs->state[i].getL();
      cout << setw(3) << bs->state[i].getS();
      cout << setw(3) << bs->state[i].getJ();
      cout << setw(3) << bs->state[i].getP() << endl;
      i++;
    }
  }
*/
  return(0);
}



/**********************************************************/
/*      Conversion from Cylindrical to Spherical          */
/*      ----------                                        */
/*      L. Bonneau, et al. Phys. Rev. C 75, 054618 (2007) */
/**********************************************************/
void HOSphericalExpansion(Basis *dbase, SPEnergy *spe, Operator *H, HFInterface *hfsp)
{
  const int    sp_nmax =  5;
  const int    sp_lmax = 15;
  const double sp_emax = 10.0;
  const double eps     = 1.0e-03;
  const double epsv2   = 1.0e-05;

  /*** Spherical Harmonics Oscillator basis */
  SphericalBasis sbase;
  sbase.setLimit(sp_nmax, sp_lmax);
  HOSphericalBasis(&sbase);

  /*** for all single-particle state in Cylindrical HO */
  hfsp->n = 0;
  for(int k=0 ; k<dbase->nv ; k++){

    if(spe->energy[k] > sp_emax) break;
    if(k >= HF_MAX_SPSTATE) break;
    if(spe->v2[k] < epsv2) break;

    int k1 = spe->subindex[k];    // index inside block
    int kb = spe->block[k];       // block index

    int omega1 = dbase->omega[kb]; // 2*Omega in deformed HO

    double piave = 0.0;
    double pitot = 0.0;

    /*** states in Spherical HO */
    int nstate = 0, lmax = 0;
    for(int i=0 ; i < sbase.nv ; i++){

      if((sbase.getM(i) != omega1) || (sbase.getJ(i) < omega1)) continue;

      double sum = OverlapFunction(i,k1,kb,dbase,&sbase,H);

      double sumsq = sum*sum;
      piave += sumsq * ( (sbase.getL(i)%2 == 0) ? 1.0 : -1.0);
      pitot += sumsq;

      /*** save state information in the HFInterface object */
      if(abs(sum) > eps){
        hfsp->state[k].set(sbase.getN(i),sbase.getL(i),sbase.getJ(i),sum);
        nstate ++;
        if(sbase.getL(i) > lmax) lmax = sbase.getL(i);
        if(nstate >= HF_MAX_EXPANSION) break;
      }
    }

    /*** average parity */
    spe->parity[k] = (pitot > 0.0) ? piave / pitot : 0.0;

    /*** save s.p. energy and v2 in the HFInterface object */
    hfsp->state[k].energy = spe->energy[k];
    hfsp->state[k].parity = (spe->parity[k] > 0.0) ? 1 : -1;
    hfsp->state[k].v2     = spe->v2[k];
    hfsp->state[k].k2     = omega1;
    hfsp->state[k].lmax   = lmax;
    hfsp->n++;
  }
  hfsp->b = 1.0/dbase->b;
}


/**********************************************************/
/*      Sum of Overlap Functions                          */
/**********************************************************/
double OverlapFunction(const int i, const int k1, const int kb,
                       Basis *dbase, SphericalBasis *sbase, Operator *H)
{
  int omega0 = sbase->getM(i);
  int nperp0 = 2*sbase->getN(i) + sbase->getL(i);

  double sum = 0.0;
  for(int k2 = 0 ; k2<dbase->ndata[kb] ; k2++){
    int k0 = dbase->index[kb] + k2;

    int nz = dbase->state[k0].getNz();
    int nr = dbase->state[k0].getNr();
    int nl = dbase->state[k0].getL();  // Lambda
    int ns = dbase->state[k0].getS();  // 2 Sigma

    int nperp2 = 2*nr + nl;
    int dnperp = nperp0 - nperp2;

    if( (nz + nperp2 - nperp0)%2 !=0 ) continue;
    if(nperp2 > nperp0) continue;
    if(nl > sbase->getL(i)) continue;

    double x = sbase->getL(i) + 0.5 + sbase->getS(i) * ns * 0.5*omega0;
    double g = sqrt(x / (2.0*sbase->getL(i) + 1.0));
    if(sbase->getS(i) == -1) g *= -ns;

    double a = OverlapFunctionA(nz,dnperp,dbase->q);
    double b = OverlapFunctionB(dnperp,nperp2,nl,sbase->getN(i),sbase->getL(i));
    x = a * b * g * H[kb].vector[k1][k2];
/*
    if(x*x > 0.01){
       cout << "    " << setw(4) << k1;
       cout << setw(4) << nz << setw(4) << nperp2;
       cout << setw(4) << nl << setw(4) << ns;
       cout << setw(12) << H[kb].vector[k1][k2];
       cout << setw(4) << sbase->getN(i) << setw(4) << sbase->getL(i);
       cout << setw(4) << 2*(sbase->getJ(i) + sbase->getL(i))-3;
       cout << setw(12) << a << setw(12) << b;
       cout << setw(12) << g << setw(12) << x << endl;
    }
*/
    sum += x;
  }
  return(sum);
}


/**********************************************************/
/*      Overlap Function Eq. (A2)                         */
/*      A_{n n'}(q)                                       */
/**********************************************************/
void HOSphericalBasis(SphericalBasis *sbase)
{
  /*** count total number of states */
  sbase->nv = 0;
  for(int l=0 ; l <=sbase->getLmax() ; l++){
    for(int s=-1 ; s<=1 ; s+=2){
      int j = 2*l + s;
      if(j < 0) continue;
      for(int m=1 ; m<=j ; m+=2){
        for(int n=0 ; n<=sbase->getNmax() ; n++){
          sbase->nv ++;
          continue;
        }
      }
    }
  }

  sbase->memalloc();

  int nc = 0;
  for(int l=0 ; l <=sbase->getLmax() ; l++){
    for(int s=-1 ; s<=1 ; s+=2){
      int j = 2*l + s;
      if(j < 0) continue;
      for(int m=1 ; m<=j ; m+=2){
        for(int n=0 ; n<=sbase->getNmax() ; n++){
          sbase->state[nc++].set(n,l,s,j,m);
        }
      }
    }
  }
/*
  for(int k2=0 ; k2<sbase->nv ; k2++){
    cout << setw(5) << k2;
    cout << setw(5) << sbase->state[k2].getN();
    cout << setw(5) << sbase->state[k2].getL();
    cout << setw(5) << sbase->state[k2].getS();
    cout << setw(5) << sbase->state[k2].getJ();
    cout << setw(5) << sbase->state[k2].getM() << endl;
  }
*/
}


/**********************************************************/
/*      Overlap Function Eq. (A2)                         */
/*      A_{n n'}(q)                                       */
/**********************************************************/
double OverlapFunctionA(int n1, int n2, double q)
{
  double a = 0.0;

  if((n1-n2)%2 != 0){ return(a); }

  if(abs(q-1.0) <= 1.0e-9){
    a = 0.0;
    if(n1 == n2) a = 1.0;
    return(a);
  }

  int m_min = max(0,(n1 - n2) / 2);
  int m_max = (n1 - n1%2) / 2;

  double x1 = (q-1.0) / (4.0 * sqrt(q));

  double sum = 0.0;
  for(int m = m_min ; m<=m_max ; m++){
    double p  = (m%2 == 0) ? 1 : -1;
    double x2 = pow(x1,2*m);
    double x3 = exp( fact[m] + fact[m+(n2-n1)/2] + fact[n1-2*m] );
    sum += p*x2/x3;
  }

  double c1 = sqrt( pow(2.0, (double)(n1-n2)) * exp(fact[n1] + fact[n2]) );
  double c2 = pow( (q-1.0) / (q+1.0), (n2-n1)/2.0 );
  double c3 = pow( 2.0 * sqrt(q) / (q+1.0), n1+0.5 );

  a = c1 * c2 * c3 * sum;

  return(a);
}


/**********************************************************/
/*      Overlap Function Eq.(A5)                          */
/*      B_{nz,nperp,Lambda;n,l}                           */
/**********************************************************/
double OverlapFunctionB(int nz, int np, int lambda, int n, int l)
{
  double b = 0.0;

  if( (np-lambda)%2 != 0 ) return(b);
  if( (nz+np) != (2*n+l) ) return(b);

  int alpha = (np + lambda)/2; // Eq. (A6)
  int beta  = (np - lambda)/2; // Eq. (A7)

  double sum = 0.0;
  for(int k=0 ; k<=l ; k++){
    int k2 = k + n - beta;
    int k1 = nz - 2*k2;

    if(k1<0 || k2<0) continue;

    double p  = (k%2 == 0) ? 1 : -1;
    double x1 = exp(fact[2*(l-k)] + fact[k+n]);
    double x2 = exp(fact[k] + fact[l-k] + fact[k1] + fact[k2]);
    sum += p * x1 / x2;
  }

  double c1 = ((alpha+n)%2 == 0) ? 1.0 : -1.0;
  double c2 = pow(2.0, (double)n);
  double c3 = exp(fact[alpha] - fact[beta]);
  double c4 = exp(fact[nz] + fact[l-lambda] + fact[n+l]) * (2*l + 1.0);
  double c5 = pow(2.0, (double)nz);
  double c6 = exp(fact[l+lambda] + fact[n] + fact[2*(n+l)+1]); 
  b = sum * c1 * c2 * sqrt(c3 * c4 / (c5 * c6));

  return(b);
}


