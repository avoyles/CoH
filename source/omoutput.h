// functions for data output in optical model calculations
// these functions are independent from output.cpp

/**************************************/
/*      omoutput.cpp                  */
/**************************************/
void    outOMP                (int, Optical *);
void    outOMPtable           (int, Optical *);
void    outCoupledState       (NuclearModel, int, LevelData *);
void    outDeformation        (NuclearModel, int, double *);
void    outLevelExcite        (int, int, double, double, LevelData *, double *);
void    outRmatrix            (int, complex<double> *, complex<double> *);
void    outSmatrix            (int, int, complex<double> *);
void    outTransmission       (int, int, double, double *);
void    outCCSmatrix          (int, int, int, Collective *, CCdata *, complex<double> *);
