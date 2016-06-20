#ifndef	KNOCK_OUT_PUT_H
#define KNOCK_OUT_PUT_H

#include <cmath>
#include "GenericQUAD.h"
#include "BetterGrid.h"
#include "LVProcess.h"


using std::cout; 

// The number of standard deviations used for computing the integration range.
// Should be a variable, I guess, but can't be bothered to keep passing it as an argument. 

/////////////////////////////////////////////////////////////////////////////////////////////////

Grid NewInitGridForDOCalls(double Strike, double Barrier, double Time, LVProcess& Proc, double Step);
Grid* InitGridForDOCalls(double Strike, double Barrier, double Time, LVProcess& Proc, unsigned long N); 
// Computes the payoffs of Euro puts for N grid points. 
// Process needed for computing the grid points, range etc.

Grid NewMakeTheNextGridForDOCalls(double Strike, double Step, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid& OldGrid);
void MakeTheNextGridForDOCalls(double Strike, double Barrier, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid* OldGrid, Grid* NewGrid); 
// Auxiliary routine needed in the backward induction 

double NewDOCall(double Strike, double Barrier, vector<double> Times, LVProcess& Proc, double Step);
double DOCall(double Strike, double Barrier, vector<double> Times, LVProcess& Proc, unsigned long N); 
// Similar routine for pricing Bermudan puts with observation points given in Times

/////////////////////////////////////////////////////////////////////////////////////////////////

double NewDOCall(double Strike, double Barrier, vector<double> Times, LVProcess& Proc, double Step)
{
	double Time = Times.back();

	Grid OldGrid = NewInitGridForDOCalls(Strike,Barrier,Time,Proc,Step);

	int i = Times.size()-1; 

	double NewTime, TimeStep, NewLow, NewHigh;

	// Work inductively backwards: while we're not at the final step, 
	// keep computing the max(Exercise,Don't)-type values 
	while (i > 0)
	{
		// Time up to the i^th observation point. Needed for computing the integration range.
		NewTime = Times[i-1]; 
		TimeStep = Times[i]-Times[i-1]; 

		// Find provisional limits for the integration
		NewLow = max(Barrier,Proc.InitLow(NewTime,DEVIATIONS)); 
		NewHigh = max(NewLow,Proc.InitHigh(NewTime,DEVIATIONS)); 
		
		// Compute the next grid store it into ReplacableGrid
		OldGrid = NewMakeTheNextGridForDOCalls(Strike,Step,NewLow,NewHigh,TimeStep,Proc,OldGrid);
		i--;
	}

	// Finally, compute the value of the option from the values at the first observation point
	double Result = SubIntegral(Proc.GetValue(),Times[0],Proc,&OldGrid); 

	return Result;
}

Grid NewInitGridForDOCalls(double Strike, double Barrier, double Time, LVProcess& Proc, double Step)
{
	// Work out the integration range and make large enough a grid
	double Low = max(Strike,Proc.InitLow(Time,DEVIATIONS)); 
	double High = max(Low,Proc.InitHigh(Time,DEVIATIONS)); 
	unsigned long N = static_cast<unsigned long>(ceil((High-Low)/Step));
	double *Values = new double[2*N+1]; 

	for (unsigned long i = 0; i < 2*N+1; i++)
	{
		Values[i] = max((Low+i*Step/2.0)-Strike,0.0); 
	}

	// Make these objects into a Grid.
	// NB: The Step size in Grid is *half* the Simpson type step size
	Grid Result = Grid(Low,Step/2.0,2*N+1,Values); 

	delete [] Values;

	return Result;
}

Grid NewMakeTheNextGridForDOCalls(double Strike, double Step, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid& OldGrid)
{
	unsigned long N = static_cast<unsigned long>(ceil((NewHigh-NewLow))/Step);
	double *Values = new double[2*N+1];

	Grid Result = Grid(NewLow,Step/2.0,2*N+1,Values); 

	for (unsigned long i = 0; i < 2*N+1; i++)
	{
		double Val = NewLow+i*Step/2.0;
		Result.SetNthValue(i,SubIntegral(Val,TimeStep,Proc,&OldGrid)); 
	}

	return Result;
}


double DOCall(double Strike, double Barrier, vector<double> Times, LVProcess& Proc, unsigned long N)
{
	// Total lifetime of the option
	double Time = Times.back(); 

	// Compute the initial grid 
	Grid* OldGrid = InitGridForDOCalls(Strike,Barrier,Time,Proc,N); 
	Grid* ReplaceableGrid = OldGrid->clone(); 
	Grid* tmp = OldGrid->clone(); 

	int i = Times.size()-1; 

	double NewTime, TimeStep, NewLow, NewHigh, NewStep;

	// Work inductively backwards: while we're not at the final step, 
	// keep computing the max(Exerise,Don't)-type values 
	while (i > 0)
	{
		// Time up to the i^th observation point. Needed for computing the integration range.
		NewTime = Times[i-1]; 
		TimeStep = Times[i]-Times[i-1]; 
		NewLow = max(Proc.InitLow(NewTime,DEVIATIONS),Barrier); 
		NewHigh = Proc.InitHigh(NewTime,DEVIATIONS); 
		NewStep = (NewHigh-NewLow)/static_cast<double>(2*N); 
		
		// The payoff is smooth above Barrier, so there's no need to adjust the grid
		// NewLow = AdjustNewLow(NewLow,Strike,NewStep); 
		// NewHigh  = NewLow + 2*N*NewStep;

		// Compute the next grid store it into ReplacableGrid
		MakeTheNextGridForDOCalls(Strike,Barrier,NewLow,NewHigh,TimeStep,Proc,OldGrid,ReplaceableGrid);

		// Swap the addresses of OldGrid & NewGrid, so that New becomes old
		Grid* tmp = OldGrid;
		OldGrid = ReplaceableGrid;
		ReplaceableGrid = tmp;

		i--;
	}

	// Finally, compute the value of the option from the values at the first observation point
	double Result = SubIntegral(Proc.GetValue(),Times[0],Proc,OldGrid); 

	// Clean up
	delete tmp;
	delete ReplaceableGrid;
	delete OldGrid;

	return Result;
}

Grid* InitGridForDOCalls(double Strike, double Barrier, double Time, LVProcess& Proc, unsigned long N)
// Creates the payoffs for a Euro put
{
	// No point in evaluating below strike, so we take...
	double Low = max(Strike,Proc.InitLow(Time,DEVIATIONS)); 
	double High = Proc.InitHigh(Time,DEVIATIONS); 
	double Step = (High-Low)/static_cast<double>(2*N); 

	double *Values = new double[2*N+1]; 

	for (unsigned long i = 0; i < 2*N+1; i++)
	{
		Values[i] = max(Low+i*Step-Strike,0.0); 
	}

	// Make these objects into a Grid
	Grid* Result = new Grid(Low,Step,2*N+1,Values); 

	delete [] Values;

	return Result;
}

void MakeTheNextGridForDOCalls(double Strike, double Barrier, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid* OldGrid, Grid* NewGrid)
{
	unsigned long N = OldGrid->GetSize(); 

	double Step = (NewHigh - NewLow)/(static_cast<double>(N) - 1.0); 
	
	NewGrid->SetLow(NewLow); 
	NewGrid->SetStep(Step);

	for (unsigned long i = 0; i < N; i++)
	{
		NewGrid->SetNthValue(i,SubIntegral(NewLow+i*Step,TimeStep,Proc,OldGrid)); 
	}
}

#endif