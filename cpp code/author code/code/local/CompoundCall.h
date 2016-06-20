#ifndef COMPOUND_CALL_H
#define COMPOUND_CALL_H

#include <cmath>
#include "GenericQUAD.h"
#include "BetterGrid.h"
#include "LVProcess.h"


using std::cout; 

// The number of standard deviations used for computing the integration range.
// Should be a variable, I guess, but can't be bothered to keep passing it as an argument. 

/////////////////////////////////////////////////////////////////////////////////////////////////

Grid InitGridForCompoundCalls(double Strike, double Time, LVProcess& Proc, unsigned long N); 
// Computes the payoffs of Euro puts for N grid points. 
// Process needed for computing the grid points, range etc.

Grid MakeTheNextGridForCompoundCalls(double Strike, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid& OldGrid); 
// Auxiliary routine needed in the backward induction 

double CompoundCall(vector<double> Strikes, vector<double> Times, LVProcess& Proc, unsigned long N); 
// Similar routine for pricing Bermudan puts with observation points given in Times

// double FloatingCompoundCall(double BaseStrike, vector<double> Times, LVProcess& Proc, unsigned long N); 
/////////////////////////////////////////////////////////////////////////////////////////////////



Grid InitGridForCompoundCalls(double Strike, double Time, LVProcess& Proc, unsigned long N)
// Creates the payoffs for a Euro put
{
	// Compute the range of the grid. The "min" appears as there's no point to 
	// add a sequence of zeros into the grid
	double Low = max(Strike,Proc.InitLow(Time,DEVIATIONS)); 
	double High = max(Strike,Proc.InitHigh(Time,DEVIATIONS));

	double Step = (High-Low)/(2.0*N+1.0); 

	double *Values = new double[2*N+1]; 
	
	for (unsigned long i = 0; i < 2*N+1; i++)
	{
		double Spot = Low+i*Step; 
		Values[i] = Spot-Strike; 
	}
		// Make these objects into a Grid
	Grid Result(Low,Step,2*N+1,Values); 
	return Result;
}

Grid MakeTheNextGridForCompoundCalls(double Strike, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid& OldGrid)
{
	unsigned long N = OldGrid.GetSize(); 

	double Step = (NewHigh - NewLow)/static_cast<double>(N); 

	double *NewGridVals = new double[N];

	for (unsigned long i = 0; i < N; i++)
	{
		// Compute the new log(Value)
		double Val = NewLow + i*Step;

		// Store the price of the remaining option in the grid
		NewGridVals[i] = max(SubIntegral(Val,TimeStep,Proc,OldGrid)-Strike,0.0);
	}

	// Make a Grid by recording the NewLow and the Step size, and return
	Grid Result(NewLow,Step,N,NewGridVals);

	return Result;
}

double CompoundCall(vector<double> Strikes, vector<double> Times, LVProcess& Proc, unsigned long N)
{
	// Total lifetime of the option
	double Time = Times.back(); 

	// Compute the initial grid 
	Grid OldGrid = InitGridForCompoundCalls(Strikes.back(),Time,Proc,N); 

	int i = Times.size()-1; 


	// Work inductively backwards: while we're not at the final step, 
	// keep computing the max(Exerise,Don't)-type values 
	while (i > 0)
	{
		// Time up to the i^th observation point. Needed for computing the integration range.
		double NewTime = Times[i-1]; 
		double NewStrike = Strikes[i-1]; 
		double TimeStep = Times[i]-Times[i-1]; 
		double NewLow = max(NewStrike,Proc.InitLow(NewTime,DEVIATIONS)); 
		double NewHigh = max(NewStrike,Proc.InitHigh(NewTime,DEVIATIONS)); 

		// Compute the next grid and replace
		Grid NewGrid = MakeTheNextGridForCompoundCalls(NewStrike,NewLow,NewHigh,TimeStep,Proc,OldGrid); 
		OldGrid = NewGrid; 
		i--;
	}

	// Finally, compute the value of the option from the values at the first observation point
	return SubIntegral(Proc.GetValue(),Times[0],Proc,OldGrid); 
}

/*
double FloatingCompoundCall(double BaseStrike, vector<double> Times, LVProcess& Proc, unsigned long N)
{
	// Total lifetime of the option
	double Time = Times.back(); 

	// Compute the initial grid 
	Grid OldGrid = InitGridForCompoundCalls(BaseStrike,Time,Proc,N); 

	int i = Times.size()-1; 


	while (i > 0)
	{
		// Time up to the i^th observation point. Needed for computing the integration range.
		double NewTime = Times[i-1]; 
		double TimeStep = Times[i]-Times[i-1]; 
		double InitVol = Proc.InitVol(); 
		double NewLow = max(NewStrike,Proc.Spot*exp((Proc.Rate-0.5*InitVol*InitVol)*NewTime - 10.0*sqrt(NewTime)*InitVol)); 
		double NewHigh = max(NewStrike,NewLow*exp(10.0*sqrt(NewTime)*InitVol)); 

		// Compute the next grid and replace
		Grid NewGrid = MakeTheNextGridForCompoundCalls(NewStrike,NewLow,NewHigh,TimeStep,Proc,OldGrid); 
		OldGrid = NewGrid; 
		i--;
	}

	// Finally, compute the value of the option from the values at the first observation point
	return SubIntegral(Proc.GetValue(),Times[0],Proc,OldGrid); 
}
*/
#endif