/*
   terminate.h : 
        prototype of function to terminate code execution
        the function is defined in main.cpp
 */

#ifndef COH_TOPLEVEL
#include <sstream>
extern ostringstream  message;
#endif

/**************************************/
/*      main.cpp                      */
/**************************************/
int     cohTerminateCode      (std::string);
void    cohNotice             (std::string);
