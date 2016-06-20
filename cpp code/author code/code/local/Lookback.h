#ifndef LOOKBACK_H
#define LOOKBACK_H

//const double DEVIATIONS = 7.5;

#include "LVProcess.h"
#include "SquareGrid.h"
#include "GenericQUAD.h"

double LookbackCall(double Strike,vector<double> Times,LVProcess& Proc, unsigned long N);

SquareGrid MakeFinalGrid(double Time, double Strike, LVProcess& Proc, unsigned long N);

SquareGrid MakeTheNextGrid(double TimeStep, double Strike, LVProcess& Proc, Grid& OldGrid); 

//----------------------------------------------------------------------------------------

SquareGrid MakeFinalGrid(double Time, double Strike, LVProcess& Proc, unsigned long N)
{
	double Spot = Proc.GetValue(); 
	double Low = Proc.InitLow(Time,DEVIATIONS);
	double High = Proc.InitHigh(Time,DEVIATIONS); 
	double Step = (High-Low)/static_cast<double>(N);

	if (Strike = Spot)
	{
		double M = ceil((Spot-Low)/Step);
		Low = Spot - M*Step; 
	}
	else
	{
		double Big = max(Strike,Spot);
		double Small = min(Strike,Spot); 
		double M = ceil((Big-Small)/Step);
		Step = (Big-Small)/M; 
		M = ceil((Small-Low)/Step);
		Low = Small - M*Step;
	}

	double **NewValues;
	NewValues = new double*[N];

	for (unsigned long i = 0; i < N; i++)
	{
		double *MVector;
		MVector = new double[N]; 
		double MValue = Low + i*Step; 
		for (unsigned long j = 0; j < N; j++)
		{
			double Value = max(MValue,Low + j*Step); 
			MVector[j] = max(Value-Strike,0.0); 
		}
		NewValues[i] = MVector; 
	}

	SquareGrid Grid(Low,Step,N,NewValues); 
	return Grid; 
}

SquareGrid MakeTheNextGrid(double TimeStep, double Strike, LVProcess& Proc, SquareGrid& OldGrid)
{
	unsigned long N = OldGrid.GetSize();
	double Low = OldGrid.GetLow();
	double Step = OldGrid.GetStep(); 

	SquareGrid NewGrid(Low,Step,N); 

	for (unsigned long i = 0; i < N; i++)
	{
		for( unsigned long j = 0; j  < N; j++)
		{
			double Param = max(Low+i*Step,Low+j*Step); 
			double NewVal = IntegrateRow(Param,TimeStep,Proc,OldGrid,i); 
			NewGrid.SetValue(i,j,NewVal); 
		}
	}

	return NewGrid; 
}

double LookbackCall(double Strike,vector<double> Times,LVProcess& Proc, unsigned long N)
{
	double FinalTime = Times.back(); 
	SquareGrid  OldGrid = MakeFinalGrid(FinalTime,Strike,Proc,N);

	unsigned long Steps = Times.size() - 1; 

	while (Steps > 0)
	{
		double TimeStep = Times[Steps] - Times[Steps-1]; 
		SquareGrid NewGrid = MakeTheNextGrid(TimeStep,Strike,Proc,OldGrid); 
		Steps = Steps - 1;
		OldGrid = NewGrid;
	}
	
	unsigned long M = static_cast<unsigned long>((Proc.GetValue() - OldGrid.GetLow())/OldGrid.GetStep());

	return IntegrateRow(Proc.GetValue(),Times[0],Proc,OldGrid,M); 
}

#endif