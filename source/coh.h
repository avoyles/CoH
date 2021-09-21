/*
   coh.h : prototype of top level functions
 */

// Define each section-head to be read in a file
enum section{BEGIN,DATA,DWBA,DSD,PREEQ,HFS,FNS,MFT,END};
typedef enum section Section;


/**************************************/
/*      coh.cpp                       */
/**************************************/
void    cohAllocateNucleus         ();
void    cohAllocateNuclearStructure();
void    cohResetAllocationPointer  ();


/**************************************/
/*      cohmain.cpp                   */
/**************************************/
void     cohMainLoop          (double, int, int, unsigned long);


/**************************************/
/*      setup.cpp                     */
/**************************************/
void    setupInitSystem       (int, Pdata *);
void    setupClearParameter   (int, System *, Pdata *, Direct *, Dcapt *, FNSpec *);
void    setupGeneralParameter (System *, Pdata *, Direct *);
int     setupStatModel        (int, System *, Pdata *);
void    setupClearCrxArray    (void);


/**************************************/
/*      dataread.cpp                  */
/**************************************/
int     readHead              (char *);
int     readSystem            (char *, System *, Pdata *, Direct *, double *, double *);
int     readDwba              (char *, Direct *);
int     readDsd               (char *, Dcapt *);
int     readStatModel         (char *, GDR *, Fission *);
int     readExciton           (char *);
int     readFns               (char *, FNSpec *);
int     readMeanfield         (char *, MFTparm *);

/**************************************/
/*      statmodel.cpp                 */
/**************************************/
int     statModel             (System *, Pdata *, Dcapt *, Direct *, FNSpec *,
                               unsigned long);

/**************************************/
/*      mftcalc.cpp                   */
/**************************************/
int     mftCalc               (MFTparm *, HFInterface *, double *);
