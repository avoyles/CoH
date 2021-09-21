/*
   fileout.h : 
        prototype of functions for result output
 */

#define filename_cross_section        "CoHCrossSection"
#define filename_particle_production  "CoHParticleProduction"
#define filename_fission              "CoHFission"
#define filename_level_excite         "CoHLevelExcite"
#define filename_angular_distribution "CoHAngularDistribution"
#define filename_legendre_coefficient "CoHLegendreCoefficient"
#define file_extension                "dat"


/**************************************/
/*      fileout.cpp                   */
/**************************************/
void    cohFileOut (int, int, double);
