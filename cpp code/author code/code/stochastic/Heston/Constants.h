#ifndef GLOBAL_CONSTANTS
#define GLOBAL_CONSTANTS

#include <cmath>
#include <complex>

typedef std::complex<double> Complex; 

const double PI = 4.0*atan(1.0);
const Complex I(0.0,1.0); 
const double OneByPI = 0.25/atan(1.0); 

#endif