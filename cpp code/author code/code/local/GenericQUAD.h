#ifndef GENERIC_QUAD_H
#define GENERIC_QUAD_H

#include "LVProcess.h"
#include "BetterGrid.h"
#include <omp.h>

const double DEVIATIONS = 6.75; 
const double vDEVIATIONS = 0.5; 
const double TOLERANCE = 0.000000001;
const double TimeTOLERANCE = 0.1;

vector<double> SplitTime(double Time, unsigned N);

double SubIntegral(double Value, double Time, LVProcess& Proc, Grid& GridPts); 


double QUADEuroCall(double Strike, double Time, LVProcess& Proc, unsigned long N); 
// Stand-alone QUAD routine for pricing European call options for local vol processes. 
// Requires that Proc has method Density.

double FindDiscontinuity(double Strike, double TimeStep,LVProcess& Proc, Grid* OldGrid);

double AdjustNewLow(double ProvLow, double Disc, double NewStep);

//------------------------------------------------------------------------------------------------


double SubIntegral(double Value, double Time, LVProcess& Proc, Grid* GridPts)
{
	double Low = GridPts->GetLow(); 
	double High = GridPts->GetHigh(); 
	double Step = GridPts->GetStep(); 

	// Compute the new integration range
	double NewLow = max(Proc.InitLow(Value,Time,DEVIATIONS),Low); 
	double NewHigh = min(Proc.InitHigh(Value,Time,DEVIATIONS),High); 

	double Result=0.0; 

	// Check whether the integration range is empty or not, then proceed.
	if (NewHigh > NewLow)
	{
		// Compute the entries of the grid corresponding to NewLow and NewHigh
		unsigned long iMinus = static_cast<unsigned long>(floor((NewLow-Low)/Step)); 
		unsigned long iPlus = static_cast<unsigned long>(floor((NewHigh-Low)/Step));

		// Must check that the range has an odd number of values (for Simpson)

		if (((iPlus-iMinus)%2))
		{
			if (iPlus < static_cast<unsigned long>(GridPts->GetSize())-1)
			{
				iPlus += 1;
			}
			else
			{
				iMinus -= 1;
			}
		}
		
		unsigned long iN = iPlus-iMinus+1;
		double Fact[iN]  Dent[iN];
		# pragma omp parallel
		{
			# pragma omp for nowait
			for (unsigned long i = 0; i <iN ; i++)
			{
				if (i == 0 || i == iN-1)
				{
					Fact[i]=1.0;	
				}
				else 
				{
					if (i%2)
					{
						Fact[i]=4.0;
					}
					else
					{
						Fact[i]=2.0;
					}
				}
			}
			
			# pragma omp for
			for (unsigned long i = 0; i <iN ; i++)
			{
				 Dent[i]=GridPts->GetNthValue(i+iMinus)*Proc.Density(Time,Low+(i+iMinus)*Step,Value);
			}

			# pragma omp for reduction(+:Result) nowait
			for (unsigned long i = 0; i <iN ; i++)
			{
				Result += Fact[i]*Dent[i];
			}
		}
	}
	
/* serier version
		// Add the first and the last term...
		Result = GridPts->GetNthValue(iMinus)*Proc.Density(Time,Low+iMinus*Step,Value);
		Result += GridPts->GetNthValue(iPlus)*Proc.Density(Time,Low+iPlus*Step,Value); 

		// ... and then the rest
		for (unsigned long i = iMinus+1; i < iPlus; i++)
		{
			// Add weights according to Simpson's rule
			if ((i-iMinus)%2)
			{
				Result += 2.0*GridPts->GetNthValue(i)*Proc.Density(Time,Low+i*Step,Value); 
			}
			else
			{
				Result += 4.0*GridPts->GetNthValue(i)*Proc.Density(Time,Low+i*Step,Value);
			}
		}
	}
*/

	// Discount by time and multiply by factors according to Simpson's rule
	Result *= exp(-Time*Proc.GetRate())*Step/3.0; 
	return Result;
}


double FindDiscontinuity(double Strike, double TimeStep,LVProcess& Proc, Grid* OldGrid)
{
	double OldX = Strike;
	double Step = -0.01;
	double NewX = OldX + Step;
	double gOldX = Strike - OldX - SubIntegral(OldX,TimeStep,Proc,OldGrid);

	double GridStep = OldGrid->GetStep(); 
	unsigned long N = 0; 

	while (fabs(Step) > TOLERANCE && N < 15)
	{
		double gNewX = Strike - NewX - SubIntegral(NewX,TimeStep,Proc,OldGrid); 
		Step = Step*gNewX/(gOldX - gNewX);
		OldX = NewX; 
		NewX = NewX - Step;
		gOldX = gNewX; 
		N++;
	}

	return NewX;
}

double AdjustNewLow(double ProvLow, double Disc, double NewStep)
{
	if (Disc < ProvLow)
	{
		return ProvLow;
	}
	else
	{
		unsigned long i = static_cast<unsigned long>(floor((Disc-ProvLow)/NewStep));

		if ((i%2))
		{
			i = i+1;
		}
		return Disc - i*NewStep;
	}
}

vector<double> SplitTime(double Time, unsigned N)
{
	double tmp = Time/static_cast<double>(N); 

	vector<double> Result;

	for (unsigned i = 1; i < N+1; i++)
	{
		Result.push_back(i*tmp); 
	}

	return Result;
}

double QUADEuroCall(double Strike, double Time, LVProcess& Proc, unsigned long N)
// Computes the price of a European call by integrating max(S_t-Strike,0) against the density
// f of the process
{
	// Compute the bounds and the Step size for the integral; 
	// we'll sum over values f(Lower+j*Step).
	// The "max" appears as there's no point in integrating over
	// the range where the option is not exercised.

	double Lower = max(Strike,Proc.InitLow(Time,DEVIATIONS)); 
	double Upper = max(Strike,Proc.InitHigh(Time,DEVIATIONS)); 
	double Step = (Upper-Lower)/(2.0*static_cast<double>(N)); 

	// Start summing the integral

	double Result = max(Lower-Strike,0.0)*Proc.Density(Time,Lower); 
	Result += max(Upper-Strike,0.0)*Proc.Density(Time,Upper); 

	// Then add weights according to Simpson's rule and sum

	unsigned long i = 1; 

	do
	{
		Result += max(Lower+i*Step-Strike,0.0)*Proc.Density(Time,Lower+i*Step)*4.0; 
		i = i+2; 
	}
	while (i < 2*N+1); 

	i = 2; 

	do 
	{
		Result += max(Lower+i*Step-Strike,0.0)*Proc.Density(Time,Lower+i*Step)*2.0; 
		i = i+2; 
	}
	while (i < 2*N-1); 

	// Further Simpson stuff
	Result *= Step/3.0;

	// Take the present value of the payoff and return
	Result *= exp(-Proc.GetRate()*Time); 

	return Result; 
}; 


#endif