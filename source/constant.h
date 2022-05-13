/*
   constant.h :
        array size and some default constants
*/


/****************************/
/*      ARRAY SIZE          */
/****************************/

#ifdef BeoH
const int MAX_ENERGY_BIN    = 3000 ;  /* maximum energy bins                   */
const int MAX_CHANNEL       =    7 ;  /* maximum decay channel                 */
const int MAX_COMPOUND      =  200 ;  /* maximum number of compound nucleus    */
const int MAX_NUCLIDE       =   50 ;  /* number of unique nuclides in chain    */
const int MAX_DIRECT        =    1 ;  /* number of direct reaction             */
const int MAX_FISS_CHANCE   =    5 ;  /* total number of multi-chance fission  */
#else
const int MAX_ENERGY_BIN    = 1000 ;  /* maximum energy bins                   */
const int MAX_CHANNEL       =    7 ;  /* maximum decay channel                 */
/* const int MAX_COMPOUND      =  900 ;  /* maximum number of compound nucleus    */
const int MAX_COMPOUND      = 2000 ;  /* maximum number of compound nucleus    */
const int MAX_NUCLIDE       =  350 ;  /* number of unique nuclides in chain    */
const int MAX_DIRECT        =   10 ;  /* number of direct reaction             */
const int MAX_FISS_CHANCE   =   20 ;  /* total number of multi-chance fission  */
#endif

const int MAX_J             =   60 ;  /* maximum J-value  = MAX_L              */
const int MAX_LEVELS        =  100 ;  /* maximum discrete levels               */
const int MAX_GAMMA_BRANCH  =   50 ;  /* maximum gamma-ray branches            */
const int MAX_MULTIPOL      =    5 ;  /* multipolarity, E1, M1, E2, M2, E3     */
const int MAX_GDR           =    8 ;  /* maximum GDRs and pygmy resonance      */
const int MAX_LAMBDA        =    5 ;  /* lambda_max = (MAX_LAMBDA-1)*2         */
const int MAX_ANGDIST       =  180 ;  /* maximum angular distribution points   */
const int MAX_ANGDISTLEVELS =   50 ;  /* maximum levels angular dist. given    */
const int MAX_HUMP          =    3 ;  /* maximum number of fission humps       */
const int MAX_KBAND         =    5 ;  /* maximum number of transition states   */
const int MAX_EINCIDENT     =  100 ;  /* maximum number of incident energies   */
const int MAX_MACSTEMP      =   50 ;  /* maximum number of MACS temperature    */
const int MAX_GAMMALINE     = 5000 ;  /* maximum number of gamma lines printed */


/****************************/
/*      DEFAULT PARAMETER   */
/****************************/

const double ENERGY_BIN     =  0.1 ; /* default energy bin in the continuum   */
const double NORM_FACT      = 10.0 ; /* conversion factor to [mb]             */
const int    MAX_PARTICLE   =   10 ; /* number of emitted particles in chain  */
const int    ANGLE_STEP     =    5 ; /* angle step for angular dist. printing */


