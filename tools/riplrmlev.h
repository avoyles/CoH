
const int MAX_LEVELS        =  1000;
const int MAX_GAMMA_BRANCH  =   100;


/**********************************************************/
/*   Z and A numbers of Nucleus                           */
/**********************************************************/
class ZAnumber{ 
 private:
    unsigned int Z;
    unsigned int A;
 public:
    ZAnumber(){
      Z = 0;
      A = 0;
    }
    ZAnumber(int z, int a){
      Z = z;
      A = a;
    }
    void setZA(int z, int a){
      Z = z;
      A = a;
    }
    unsigned int getZ(){ return (Z); }
    unsigned int getA(){ return (A); }
    unsigned int getN(){ return (A-Z); }
    ZAnumber operator+(ZAnumber x){
      ZAnumber y;
      y.Z = Z + x.Z;
      y.A = A + x.A;
      return y;
    }
    ZAnumber operator-(ZAnumber x){
      ZAnumber y;
      y.Z = Z - x.Z;
      y.A = A - x.A;
      return y;
    }
    bool operator==(ZAnumber x){
      bool z = false;
      if( (Z == x.Z) && (A == x.A) ) z = true;
      return z;
    }
};


/**********************************************************/
/*    Discrete Levels, energy, spin, parity               */
/**********************************************************/
class Level{
 public:
    int        flag          ;     // special flag not defined in RIPL
    double     energy        ;     // level energy
    double     spin          ;     // level spin
    double     halflife      ;     // half-life of the state
    int        parity        ;     // parity
    int        ngamma        ;     // number of Gamma-rays
    int       *fs            ;     // final state index for g-decay
    double    *eg            ;     // gamma-ray energy
    double    *pg            ;     // probability of gamma decay
    double    *pe            ;     // probability of electromagnetic transition
    double    *cc            ;     // internal conversion coefficients
    string     info          ;

    Level(){
      flag     = 0;
      energy   = 0.0;
      spin     = 0.0;
      halflife = 0.0;
      parity   = 0;
      ngamma   = 0;
      info     = "";
    }
};


/**********************************************************/
/*    Nuclear Structure Data, Discrete Levels             */
/**********************************************************/
class RIPLLevel{
 private:
    int        nlevmax       ;     // max number of levels allocated
    bool       arrayalloc    ;     // flag for memory allocation
 public:
    string     name          ;
    ZAnumber   za            ;     // Z and A numbers
    Level     *lev           ;     // discrete level data
    int        nol           ;     // number discrete levels
    int        nog           ;     // number gamma lines
    int        nmax          ;
    int        nc            ;
    double     sn            ;
    double     sp            ;

    RIPLLevel(){
      name       = "     ";
      za.setZA(0,0);
      nlevmax    = 0;
      nol        = 0;
      nog        = 0;
      nmax       = 0;
      nc         = 0;
      sn         = 0.0;
      sp         = 0.0;
      arrayalloc = false;
    }
    void memalloc(int n, int mb){
      if(!arrayalloc){
        nlevmax = n;
        lev = new Level [nlevmax];
        for(int j=0 ; j<nlevmax ; j++){
          lev[j].fs = new int[mb];
          lev[j].eg = new double[mb];
          lev[j].pg = new double[mb];
          lev[j].pe = new double[mb];
          lev[j].cc = new double[mb];
        }
        arrayalloc = true;
      }
    }
    ~RIPLLevel(){
      if(arrayalloc){
        for(int j=0 ; j<nlevmax ; j++){
          delete [] lev[j].fs;
          delete [] lev[j].eg;
          delete [] lev[j].pg;
          delete [] lev[j].pe;
          delete [] lev[j].cc;
        }
        delete [] lev;
        arrayalloc = false;
      }
    }
};



void RIPLRemoveLevel(RIPLLevel *);
bool RIPLCheckRemoveFlag(RIPLLevel *);
void RIPLPrintLevel(RIPLLevel *);
