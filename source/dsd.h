// functions for DSD model, and constants

#define MAX_POINTS_BW   400  /* maximum radial points for bound wave   */
#define MAX_SPLEVEL      50  /* maximum number of sp states            */
#define VIB_FORMFACTOR    2  /* index of coupling formfactor           */


/****************************/
/*   Bound State Data       */
/****************************/
struct boundstate{
  public:
    int     n                ;     /* bound state principle q.num. n   */
    int     l                ;     /* bound state q.num. l             */
    int     j2               ;     /* bound state q.num. j*2           */
    int     k2               ;     /* band K*2 (or Omega)              */
    double  bind             ;     /* bound state binding energy       */
    double  w                ;     /* BCS theory U factor or weight    */
    double  *wave            ;     /* bound wave func.                 */
};
typedef struct boundstate Bstate;

struct phstate{
  public:
    Bstate  Z                ;     /* proton shell                     */
    Bstate  N                ;     /* neutron shell                    */
};
typedef struct phstate PHstate;


/**************************************/
/*      BCS.CPP                       */
/**************************************/
double  bcsEnergyGap          (Particle, int, int);
double  bcsVfactor            (double, double, double);
double  bcsGapEquation        (int, double, double, Bstate *);


/**************************************/
/*      NILSSON.CPP                   */
/**************************************/
int     nilssonFillOrbit      (Particle, int);
int     nilssonLevel          (Particle, int, int, unsigned int, Bstate *, double, double, double);


/**************************************/
/*      BSTATWAV.CPP                  */
/**************************************/
int     bstateWavefunc        (Pdata *, ZAnumber *, double, double, double, double, double, 
                               Optical *, Potential *, Bstate *);


/**************************************/
/*      BSTATSOLVE.CPP                */
/**************************************/
double  bstateFindDepth       (int, int, int, double, double, double, double *,
                               Potential *);


/**************************************/
/*      DSDOVERLAP.C                  */
/**************************************/
void    dsdFormfactor         (int, Optical *, Potential *, double, double,
                               complex<double> *, complex<double> *);
complex<double> dsdRadialIntegral     (int, double, complex<double> *, double *, complex<double> *);
