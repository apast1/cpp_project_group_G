#ifndef MISC_MATHS_H
#define MISC_MATHS_H

// Miscellaneous mathematical functions apparently not in the standard library

#include "MersenneTwister.h"
#include <vector>

typedef std::complex<double> Complex; 

using std::vector; 

double max(double a, double b)
{
	if (a >= b )
	{
		return a;
	}
	else 
	{
		return b; 
	}
}

double min(double a, double b)
{
	if (a >=b)
	{
		return b;
	}
	else
	{
		return a;
	}
}

double Diff(double a, double b)
{
	return max(a-b,b-a);
}

double Normal(const double x)
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

unsigned long Power(unsigned long N, int n)
{
	unsigned Result = 1;
	unsigned long Lim = static_cast<unsigned long>(n); 

	if (n == 0)
	{
		return Result;
	}
	else 
	{
		Result = N; 

		for (unsigned long i=0; i < Lim; i++)
		{
			Result *= N; 
		}
	}

	return Result; 
}

std::vector<double> GetCorrelatedNormals(double Corr, MTRand& MT)
{
	vector<double> Result;
	Result.push_back(MT.randNorm()); 
	Result.push_back(sqrt(1-Corr*Corr)*MT.randNorm()+Corr*Result[0]); 
	return Result;
}

std::vector<double> realpart(std::vector<Complex> Cxvec) 
{
	std::vector<double> Result; 
	for (unsigned i=0; i < Cxvec.size(); i++)
	{
		Result.push_back(real(Cxvec[i]));
	}
	return Result; 
}

void Multiply(std::vector<double>& Vec, double x)
// Multiply a vector Vec by a scalar x
{
	for (unsigned long i = 0; i < Vec.size(); i++)
	{
		Vec[i] *= x; 
	}
}
#endif 