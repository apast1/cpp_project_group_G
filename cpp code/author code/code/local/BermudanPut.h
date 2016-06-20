#ifndef	BERMUDAN_PUT_H
#define BERMUDAN_PUT_H

#include <cmath>
#include "GenericQUAD.h"
#include "BetterGrid.h"
#include "LVProcess.h"


using std::cout; 

// The number of standard deviations used for computing the integration range.
// Should be a variable, I guess, but can't be bothered to keep passing it as an argument. 

/////////////////////////////////////////////////////////////////////////////////////////////////

double NewBermudanPut(double Strike, vector<double> Times, LVProcess& Proc, double Step);
// Bermudan puts with a given Strike and observation points in Times.
// NB: The step size is fixed at all observation levels; "N" varies.

double BermudanPut(double Strike, vector<double> Times, LVProcess& Proc, unsigned long N); 
// Bermudan puts with observation points given in Times. 
// NB: Uses the same number N of grid points at each level!

Grid* InitGridForPuts(double Strike, double Time, LVProcess& Proc, unsigned long N); 
// Computes the payoffs of Euro puts for N grid points. 
// Process needed for computing the grid points, range etc.


void MakeTheNextGridForPuts(double Strike, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid* OldGrid, Grid* NewGrid); 
// Auxiliary routine needed in the backward induction 

Grid NewInitGridForPuts(double Strike, double Time, LVProcess& Proc, double Step);

Grid NewMakeTheNextGridForPuts(double Strike, double Step, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid& OldGrid);


/////////////////////////////////////////////////////////////////////////////////////////////////

double NewBermudanPut(double Strike, vector<double> Times, LVProcess& Proc, double Step)
{
	double Time = Times.back();

	Grid OldGrid = NewInitGridForPuts(Strike,Time,Proc,Step);

	int i = Times.size()-1; 

	double NewTime, TimeStep, NewLow, NewHigh, Disc;

	// Work inductively backwards: while we're not at the final step, 
	// keep computing the max(Exercise,Don't)-type values 
	while (i > 0)
	{
		// Time up to the i^th observation point. Needed for computing the integration range.
		NewTime = Times[i-1]; 
		TimeStep = Times[i]-Times[i-1]; 

		// Find provisional limits for the integration
		NewLow = Proc.InitLow(NewTime,DEVIATIONS); 
		NewHigh = Proc.InitHigh(NewTime,DEVIATIONS); 
		
		// Find the point where the payoff isn't smooth, then adjust the grid accordingly
		Disc = FindDiscontinuity(Strike,TimeStep,Proc,&OldGrid);
		double Error = NewLow - AdjustNewLow(NewLow,Disc,Step/2.0); 
		NewLow = NewLow - Error;
		NewHigh  = NewHigh - Error;

		// Compute the next grid store it into ReplacableGrid
		OldGrid = NewMakeTheNextGridForPuts(Strike,Step,NewLow,NewHigh,TimeStep,Proc,OldGrid);
		i--;
	}

	// Finally, compute the value of the option from the values at the first observation point
	double Result = SubIntegral(Proc.GetValue(),Times[0],Proc,&OldGrid); 

	return Result;
}

Grid NewInitGridForPuts(double Strike, double Time, LVProcess& Proc, double Step)
{
	// Work out the integration range and make large enough a grid
	double Low = Proc.InitLow(Time,DEVIATIONS); 
	double High = min(Strike,Proc.InitHigh(Time,DEVIATIONS)); 
	unsigned long N = static_cast<unsigned long>(ceil((High-Low)/Step));
	double *Values = new double[2*N+1]; 

	// Adjust the step so that the non-linearity point is still hit
	double NewStep = (High-Low)/static_cast<double>(N); 

	for (unsigned long i = 0; i < 2*N+1; i++)
	{
		Values[i] = max(Strike-(Low+i*NewStep/2.0),0.0); 
	}

	// Make these objects into a Grid.
	// NB: The Step size in Grid is *half* the Simpson type step size
	Grid Result = Grid(Low,NewStep/2.0,2*N+1,Values); 

	delete [] Values;

	return Result;
}

Grid NewMakeTheNextGridForPuts(double Strike, double Step, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid& OldGrid)
{
	unsigned long N = static_cast<unsigned long>(ceil((NewHigh-NewLow))/Step);
	double *Values = new double[2*N+1];

	Grid Result = Grid(NewLow,Step/2.0,2*N+1,Values); 

	for (unsigned long i = 0; i < 2*N+1; i++)
	{
		// Compute the new asset value
		double Val = NewLow + i*Step/2.0;

		// Compute the value for exercising/not exercising the option
		double DontExercise = SubIntegral(Val,TimeStep,Proc,&OldGrid);

		// Add the maximum to the vector
		Result.SetNthValue(i,max(Strike-Val,DontExercise)); 
	}

	return Result;
}


double BermudanPut(double Strike, vector<double> Times, LVProcess& Proc, unsigned long N)
{
	// Total lifetime of the option
	double Time = Times.back(); 

	// Compute the initial grid 
	Grid* OldGrid = InitGridForPuts(Strike,Time,Proc,N); 
	Grid* ReplaceableGrid = OldGrid->clone(); 
	Grid* tmp = OldGrid->clone(); 

	int i = Times.size()-1; 

	double NewTime, TimeStep, NewLow, NewHigh, NewStep, Disc;

	// Work inductively backwards: while we're not at the final step, 
	// keep computing the max(Exerise,Don't)-type values 
	while (i > 0)
	{
		// Time up to the i^th observation point. Needed for computing the integration range.
		NewTime = Times[i-1]; 
		TimeStep = Times[i]-Times[i-1]; 
		NewLow = Proc.InitLow(NewTime,DEVIATIONS); 
		NewHigh = Proc.InitHigh(NewTime,DEVIATIONS); 
		NewStep = (NewHigh-NewLow)/static_cast<double>(2*N); 
		
		// Find the point where the payoff isn't smooth
		Disc = FindDiscontinuity(Strike,TimeStep,Proc,OldGrid);
		NewLow = AdjustNewLow(NewLow,Disc,NewStep); 
		NewHigh  = NewLow + 2*N*NewStep;

		// Compute the next grid store it into ReplacableGrid
		MakeTheNextGridForPuts(Strike,NewLow,NewHigh,TimeStep,Proc,OldGrid,ReplaceableGrid);

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

Grid* InitGridForPuts(double Strike, double Time, LVProcess& Proc, unsigned long N)
// Creates the payoffs for a Euro put
{
	// Compute the range of the grid. The "min" appears as there's no point to 
	//0 add a sequence of zeros into the grid
	double Low = Proc.InitLow(Time,DEVIATIONS); 
	double High = min(Strike,Proc.InitHigh(Time,DEVIATIONS)); 
	double Step = (High-Low)/static_cast<double>(2*N); 

	double *Values = new double[2*N+1]; 

	for (unsigned long i = 0; i < 2*N+1; i++)
	{
		Values[i] = max(Strike-(Low+i*Step),0.0); 
	}

	// Make these objects into a Grid
	Grid* Result = new Grid(Low,Step,2*N+1,Values); 

	delete [] Values;

	return Result;
}

void MakeTheNextGridForPuts(double Strike, double NewLow, double NewHigh, double TimeStep, LVProcess& Proc, Grid* OldGrid, Grid* NewGrid)
{
	unsigned long N = OldGrid->GetSize(); 

	double Step = (NewHigh - NewLow)/(static_cast<double>(N) - 1.0); 
	
	NewGrid->SetLow(NewLow); 
	NewGrid->SetStep(Step);

	for (unsigned long i = 0; i < N; i++)
	{
		// Compute the new log(Value)
		double Val = NewLow + i*Step;

		// Compute the value for exercising/not exercising the option
		double DontExercise = SubIntegral(Val,TimeStep,Proc,OldGrid);

		// Add the maximum to the vector
		NewGrid->SetNthValue(i,max(Strike-Val,DontExercise)); 
	}
}


#endif