/*
   levden.h : 
        prototype of functions for level density formulas
 */


/**************************************/
/*      levden.cpp                    */
/**************************************/

double  ldFermiGas                      (double, double, double);
double  ldConstantTemperature           (double, double, double);

double  ldSpinCutoff                    (double, double, double, double);
double  ldShellCorrection               (double, double, double, double);

double  ldLevelDensity                  (double, double, LevelDensity *);
double  ldDensityParameter              (double, double, LevelDensity *);
double  ldSpinDistribution              (double, double, double, double);
double  ldParityDistribution            (int);
double  ldParityDistributionEdepend     (int, double, double);
double  ldLevelSpinCutoff               (int, Level *);
void    ldTGconnect                     (double, double, double, LevelDensity *);
void    ldTextrapolate                  (double, double, LevelDensity *);
void    ldGextrapolate                  (double, LevelDensity *);
