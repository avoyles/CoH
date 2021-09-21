/*****************************************************/
/*         User defined OMPs                         */
/*****************************************************/

#include <string>
using namespace std;
#include "omplib.h"

unsigned int OMPscratch(double e,int at,int zt,int ai,int zi,Optical *omp)
{
  unsigned int pf = OMPKoningDelaroche(e,at,zt,ai,zi,omp);

  omp->a0 *= 0.6;

  return(pf);
}
