#ifndef	QUAD_H
#define QUAD_H

#include <cmath>
#include "Grid.h"

const double Deviations = 7.5; 
// The number of standard deviations used for computing the integration range.
// Should be a variable, I guess, but can't be bothered to keep passing it as an argument. 

using std::vector; 

/////////////////////////////////////////////////////////////////////////////////////////////////

Grid InitGridForPuts(double Strike, double Time, Process& Proc, unsigned long N); 
// Computes the payoffs of Euro puts for N grid points. 
// Process needed for computing the grid points, range etc.

Grid MakeTheNextGridForPuts(double Strike, double NewLow, double NewHigh, double TimeStep, Process& Proc, Grid& OldGrid); 
// Auxiliary routine needed in the backward induction 

double SubIntegral(double logVal, double Time, Process& Proc, Grid& GridPts); 
// Auxiliary function for computing integrals needed for grid points.
// More precisely, integrates over a suitable range of values given in GridPts.

double QUADEuroCall(double Strike, double Time, Process& Proc, unsigned long N); 
// Stand-alone QUAD routine for pricing European call options. 
// Requires that Proc has method Density (implemented)
// N = Number of grid points; 30+ recommended, but depends on the process

double BermudanPut(double Strike, vector<double> Times, Process& Proc, unsigned long N); 
// Similar routine for pricing Bermudan puts with observation points given in Times

vector<double> SplitTime(double Time, unsigned N); 
// Auxiliary algorithm for producing evenly spaced time grids (so that can approximate
// the price of American options by Bermudans). 

/////////////////////////////////////////////////////////////////////////////////////////////////



double QUADEuroCall(double Strike, double Time, Process& Proc, unsigned long N)
// Computes the price of a European call by integrating max(S_t-Strike,0) against the density
// f of the process
{
	double logStr = log(Strike); 
	double logS = log(Proc.GetValue()); 
	double Rate = Proc.GetRate(); 
	double Vol = Proc.GetVol(); 

	// Compute the bounds and the Step size for the integral; 
	// we'll sum over values f(Lower+j*Step).
	// The "min" appears as there's no point in integrating over
	// the range where the option is not exercised.

	double Step = (2.0*Vol*Deviations)/static_cast<double>(2*N+1);
	double Lower =  logStr; 
	double Upper =  min(logStr,logS+Deviations*sqrt(Time)*Vol); 

	// Start summing the integral

	double Result = max(exp(Lower)-Strike,0.0)*Proc.Density(Time,Lower); 
	Result += max(exp(Upper)-Strike,0.0)*Proc.Density(Time,Upper); 

	// Then add weights according to Simpson's rule and sum

	unsigned long i = 1; 

	do
	{
		Result += max(exp(Lower+i*Step)-Strike,0.0)*Proc.Density(Time,Lower+i*Step)*4.0; 
		i = i+2; 
	}
	while (i < 2*N+1); 

	i = 2; 

	do 
	{
		Result += max(exp(Lower+i*Step)-Strike,0.0)*Proc.Density(Time,Lower+i*Step)*2.0; 
		i = i+2; 
	}
	while (i < 2*N-2); 

	// Further Simpson stuff
	Result *= Step/3.0;

	// Take the present value of the payoff and return
	Result *= exp(-Rate*Time); 

	return Result; 
}; 

Grid InitGridForPuts(double Strike, double Time, Process& Proc, unsigned long N)
// Creates the payoffs for a Euro put
{
	double logValue = log(Proc.GetValue()); 
	double Rate = Proc.GetRate();
	double Vol = Proc.GetVol(); 

	// Compute the rage of the grid. The "min" appears as there's no point to 
	// add a sequence of zeros into the grid
	double Low = logValue+(Rate-0.5*Vol*Vol)*Time-Deviations*sqrt(Time)*Vol; 
	double High = min(log(Strike),Low+2*Deviations*sqrt(Time)*Vol); 
	double Step = (High-Low)/(2.0*N+1.0); 

	vector<double> Values; 

	for (unsigned long i = 0; i < 2*N+1; i++)
	{
		double Spot = exp(Low+i*Step); 
		Values.push_back(Strike-Spot); 
	}

	// Make these objects into a Grid
	Grid Result(Low,Step,Values); 

	return Result;
}

Grid MakeTheNextGridForPuts(double Strike, double NewLow, double NewHigh, double TimeStep, Process& Proc, Grid& OldGrid)
{
	unsigned long N = OldGrid.GetSize(); 

	double Step = (NewHigh - NewLow)/static_cast<double>(N); 

	vector<double> NewGrid(N,0);

	for (unsigned long i = 0; i < N; i++)
	{
		// Compute the new log(Value)
		double LogVal = NewLow + i*Step;

		// Compute the value for exercising/not exercising the option
		double Exercise = max(Strike-exp(LogVal),0.0); 
		double Dont = SubIntegral(LogVal,TimeStep,Proc,OldGrid);

		// Add the maximum to the vector
		NewGrid[i] = max(Exercise,Dont); 
	}

	// Make a Grid by recording the NewLow and the Step size, and return
	Grid Result(NewLow,Step,NewGrid);

	return Result;
}


double BermudanPut(double Strike, vector<double> Times, Process& Proc, unsigned long N)
{
	// Total lifetime of the option
	double Time = Times.back(); 

	double logValue = log(Proc.GetValue()); 
	double Rate = Proc.GetRate();
	double Vol = Proc.GetVol();

	// Compute the initial grid 
	Grid OldGrid = InitGridForPuts(Strike,Time,Proc,N); 

	int i = Times.size()-1; 

	// Work inductively backwards: while we're not at the final step, 
	// keep computing the max(Exerise,Don't)-type values 
	while (i > 0)
	{
		// Time up to the i^th observation point. Needed for computing the integration range.
		double NewTime = Times[i-1]; 
		double TimeStep = Times[i]-Times[i-1]; 

		double NewLow = logValue + (Rate-0.5*Vol*Vol)*NewTime - Deviations*sqrt(NewTime)*Vol; 
		double NewHigh = NewLow + 2*Deviations*sqrt(NewTime)*Vol; 

		// Compute the next grid and replace
		Grid NewGrid = MakeTheNextGridForPuts(Strike,NewLow,NewHigh,TimeStep,Proc,OldGrid); 
		OldGrid = NewGrid; 
		i--;
	}

	// Finally, compute the value of the option from the values at the first observation point

	return SubIntegral(logValue,Times[0],Proc,OldGrid); 
}

double SubIntegral(double logVal, double Time, Process& Proc, Grid& GridPts)
{
	double Rate = Proc.GetRate(); 
	double Vol = Proc.GetVol(); 

	// Precompute the standard term... 
	double Stupid = (Rate-0.5*Vol*Vol)*Time; 

	double Low = GridPts.GetLow(); 
	double High = GridPts.GetHigh(); 
	double Step = GridPts.GetStep(); 

	// Compute the new integration range
	double NewLow = max(logVal+Stupid-Deviations*Vol*sqrt(Time),Low); 
	double NewHigh = min(logVal+Stupid+Deviations*Vol*sqrt(Time),High); 

	double Result=0.0; 

	// Check whether the integration range is empty or not, then proceed.
	if (NewHigh > NewLow)
	{
		// Compute the entries of the grid corresponding to NewLow and NewHigh
		int iMinus = static_cast<int>(floor((NewLow-Low)/Step)); 
		int iPlus = static_cast<int>(floor((NewHigh-Low)/Step));

		// Add the first term...
		Result = GridPts.GetNthValue(iMinus)*Proc.NewDensity(Time,Low+iMinus*Step,logVal); 

		// ... and then the rest
		for (int i = iMinus+1; i < iPlus; i++)
		{
			// Add weights according to Simpson's rule
			double Factor;
			if (i%2)
			{
				Factor = 4.0;
			}
			else
			{
				Factor = 2.0;
			}
			Result += Factor*GridPts.GetNthValue(i)*Proc.NewDensity(Time,Low+i*Step,logVal); 
		}
	}

	// Discount by time and multiply by factors according to Simpson's rule
	Result *= exp(-Time*Rate)*Step/3.0; 

	return Result;
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


#endif