#ifndef SVQUAD_H
#define SVQUAD_H

#include <fstream>
#include <cmath>

using std::ofstream;

#include "GenericQUAD4.h"
#include "SquareGrid2.h"
#include "Constants.h"

double BermudanPut(double Strike, vector<double> Times, SVProcess& Proc, unsigned long Steps);

double DOCall(double Strike, double Barrier, vector<double> Times, SVProcess& Proc, unsigned long Steps);

double EuroCall(double Strike, double Time, SVProcess& Proc, unsigned long Steps);

double EuroPut(double Strike, double Time, SVProcess& Proc, unsigned long Steps);

SquareGrid InitGridForPuts(double Strike, double Time, SVProcess& Proc, unsigned long Steps); 

SquareGrid InitGridForDOCalls(double Strike, double Barrier, double Time, SVProcess& Proc, unsigned long Steps); 

SquareGrid MakeTheNextGridForPuts(double Strike,double Time,double TimeStep,SVProcess& Proc,SquareGrid& OldGrid,unsigned long Steps);

SquareGrid MakeTheNextGridForDOCalls(double Strike,double Barrier,double Time,double TimeStep,SVProcess& Proc,SquareGrid& OldGrid,unsigned long Steps);

SquareGrid MakeDummyGrid(double Time,double TimeStep,SVProcess& Proc,SquareGrid& OldGrid,unsigned long Steps);


/////////////////////////////////////////////////////////////////////////////////////////////////////////

double BermudanPut(double Strike, vector<double> Times, SVProcess& Proc, unsigned long Steps)
{
	double Time = Times.back(); 
	SquareGrid BGrid = InitGridForPuts(Strike,Time,Proc,Steps);

	int i = Times.size()-1; 
	double TimeStep = Times.front();

	while (i > 0)
	{
		double NewTime = Times[i-1];  
		BGrid = MakeTheNextGridForPuts(Strike,NewTime,TimeStep,Proc,BGrid,Steps);
		i--;
	}

//	char* grid= "gdtstPUT.csv";
//	BGrid.PrintIntoFile(grid);
	return exp(-Proc.GetRate()*Times[0])*SubIntegral(Proc.GetValue(),Proc.GetVar(),Times[0],Proc,BGrid);
}

double DOCall(double Strike, double Barrier, vector<double> Times, SVProcess& Proc, unsigned long Steps)
{
	double Time = Times.back(); 
	SquareGrid BGrid = InitGridForDOCalls(Strike,Barrier,Time,Proc,Steps);

	int i = Times.size()-1; 
	double TimeStep = Times.front();
	
	while (i > 0)
	{
		double NewTime = Times[i-1]; 
		BGrid = MakeTheNextGridForDOCalls(Strike,Barrier,NewTime,TimeStep,Proc,BGrid,Steps);
		i--;
	}


	return exp(-Proc.GetRate()*Times[0])*SubIntegral(Proc.GetValue(),Proc.GetVar(),Times[0],Proc,BGrid);
}

double EuroCall(double Strike, double Time, SVProcess& Proc, unsigned long Steps)
{
	SquareGrid BGrid = InitGridForDOCalls(Strike,0.00001,Time,Proc,Steps);
	char* grid= "gdtst.csv";
	BGrid.PrintIntoFile(grid);
	unsigned long N = static_cast<unsigned long>(ceil(Time/TimeTOLERANCE));
	double TimeStep = Time/static_cast<double>(N);

	for (unsigned long i = 1; i < N; i++)
	{
		BGrid = MakeDummyGrid(Time-i*TimeStep,TimeStep,Proc,BGrid,Steps);
	}

	return exp(-TimeStep*Proc.GetRate())*SubIntegral(Proc.GetValue(),Proc.GetVar(),TimeStep,Proc,BGrid);
}

double EuroPut(double Strike, double Time, SVProcess& Proc, unsigned long Steps)
{
	SquareGrid BGrid = InitGridForPuts(Strike,Time,Proc,Steps);
	return exp(-Time*Proc.GetRate())*SubIntegral(Proc.GetValue(),Proc.GetVar(),Time,Proc,BGrid);
}

SquareGrid InitGridForPuts(double Strike, double Time, SVProcess& Proc, unsigned long Steps)
{
	double sLow = Proc.InitsLow(Time,DEVIATIONS); 
	double sHigh = min(log(Strike),Proc.InitsHigh(Time,DEVIATIONS)); 
	double sStep = (sHigh-sLow)/static_cast<double>(Steps);

	double vLow = Proc.InitvLow(Time,vDEVIATIONS); 
	double vHigh = Proc.InitvHigh(Time,vDEVIATIONS); 
	double vStep = (vHigh-vLow)/static_cast<double>(Steps);

	SquareGrid Result(sLow,vLow,sStep,vStep,Steps,Steps);

	for (unsigned long i = 0; i < Steps+1; i++)
	{
		for (unsigned long j = 0; j < Steps+1; j++)
		{
			Result.SetValue(i,j,max(Strike-exp(sLow+i*sStep),0.0));
		}
	}

	return Result;
}

SquareGrid InitGridForDOCalls(double Strike, double Barrier, double Time, SVProcess& Proc, unsigned long Steps)
{
	double sLow = max(log(Strike),Proc.InitsLow(Time,DEVIATIONS)); 
	double sHigh = Proc.InitsHigh(Time,DEVIATIONS); 
	double sStep = (sHigh-sLow)/static_cast<double>(Steps);

	double vLow = Proc.InitvLow(Time,vDEVIATIONS); 
	double vHigh = Proc.InitvHigh(Time,vDEVIATIONS); 
	double vStep = (vHigh-vLow)/static_cast<double>(Steps);
	
	SquareGrid Result(sLow,vLow,sStep,vStep,Steps,Steps);

	for (unsigned long i = 0; i < Steps+1; i++)
	{
		for (unsigned long j = 0; j < Steps+1; j++)
		{
			Result.SetValue(i,j,max(exp(sLow+i*sStep)-Strike,0.0));
		}
	}

	return Result;
}

SquareGrid MakeTheNextGridForPuts(double Strike,double Time,double TimeStep,SVProcess& Proc,SquareGrid& OldGrid,unsigned long Steps)
{
	// Work out whether there need to be intermediate grids. Then adjust the timestep.
	// TimeTOLERANCE in a constant written down in GenericQUAD.h
	unsigned long Levels = static_cast<unsigned long>(ceil(TimeStep/TimeTOLERANCE));
	double NewTimeStep = TimeStep/static_cast<double>(Levels);

	// Set up a temp grid
	SquareGrid Old = OldGrid;

	for (unsigned long i = 1; i < Levels; i++)
	// NB: If TimeStep =< TimeTOLERANCE, this loop is empty
	{	
		// This will move back in small steps. After the loop, NewTime = Time + NewTimeStep.
		double NewTime = Time + TimeStep - i*NewTimeStep;
		Old = MakeDummyGrid(NewTime,NewTimeStep,Proc,Old,Steps);
	}
	
	double sHigh = min(log(Strike),Proc.InitsHigh(Time,DEVIATIONS));
	double sLow = Proc.InitsLow(Time,DEVIATIONS);
	double sStep = (sHigh-sLow)/static_cast<double>(Steps);

	double vHigh = Proc.InitvHigh(Time,vDEVIATIONS);
	double vLow = Proc.InitvLow(Time,vDEVIATIONS);
	double vStep = (vHigh-vLow)/static_cast<double>(Steps);

	// Make an empty grid
	SquareGrid Result(sLow,vLow,sStep,vStep,Steps,Steps);

	// Fill in the values
	for (unsigned long i = 0; i < Steps+1; i++)
	{
		for (unsigned long j = 0; j < Steps+1; j++)
		{
			double Integral = SubIntegral(sLow+i*sStep,vLow+j*vStep,NewTimeStep,Proc,Old);
			Result.SetValue(i,j,max(Strike-exp(sLow+i*sStep),Integral));
		}
	}

	
	return Result;
}

SquareGrid MakeDummyGrid(double Time,double TimeStep,SVProcess& Proc,SquareGrid& OldGrid,unsigned long Steps)
{
	double sHigh = Proc.InitsHigh(Time,DEVIATIONS);
	double sLow = Proc.InitsLow(Time,DEVIATIONS);
	double sStep = (sHigh-sLow)/static_cast<double>(Steps);		

	double vHigh = Proc.InitvHigh(Time,vDEVIATIONS);
	double vLow = Proc.InitvLow(Time,vDEVIATIONS);
	double vStep = (vHigh-vLow)/static_cast<double>(Steps);

	SquareGrid Result(sLow,vLow,sStep,vStep,Steps,Steps);

	for (unsigned long j = 0; j < Steps+1; j++)
	{
		for (unsigned long k = 0; k < Steps+1; k++)
		{
			Result.SetValue(j,k,SubIntegral(sLow+j*sStep,vLow+k*vStep,TimeStep,Proc,OldGrid));
		}
	}

	Result.PrintIntoFile("Taulukko.HJH");

	return Result;
}


SquareGrid MakeTheNextGridForDOCalls(double Strike,double Barrier,double Time,double TimeStep,SVProcess& Proc,SquareGrid& OldGrid,unsigned long Steps)
{
	// Work out whether there need to be intermediate grids. Then adjust the timestep.
	// TimeTOLERANCE in a constant written down in GenericQUAD.h
	unsigned long Levels = static_cast<unsigned long>(ceil(TimeStep/TimeTOLERANCE));
	double NewTimeStep = TimeStep/static_cast<double>(Levels);

	// Set up a temp grid
	SquareGrid Old = OldGrid;

	for (unsigned long i = 1; i < Levels; i++)
	// NB: If TimeStep =< TimeTOLERANCE, this loop is empty
	{	
		// This will move back in small steps. After the loop, NewTime = Time + NewTimeStep.
		double NewTime = Time + TimeStep - i*NewTimeStep;
		Old = MakeDummyGrid(NewTime,NewTimeStep,Proc,Old,Steps);
	}
	
	double sHigh = Proc.InitsHigh(Time,DEVIATIONS);
	double sLow = max(max(log(Strike),Proc.InitsLow(Time,DEVIATIONS)),log(Barrier));
	unsigned long NewSteps = 0;
	if (max(log(Strike),Proc.InitsLow(Time,DEVIATIONS)) < log(Barrier))
	{
		NewSteps = static_cast<unsigned long>(ceil(Steps*(sHigh-sLow)/(sHigh-max(log(Strike),Proc.InitsLow(Time,DEVIATIONS)))));
	}
	else
	{
		NewSteps = Steps;
	}
	double sStep = (sHigh-sLow)/static_cast<double>(NewSteps);

	double vHigh = Proc.InitvHigh(Time,vDEVIATIONS);
	double vLow = Proc.InitvLow(Time,vDEVIATIONS);
	double vStep = (vHigh-vLow)/static_cast<double>(Steps); 

	// Make an empty grid
	SquareGrid Result(sLow,vLow,sStep,vStep,NewSteps,Steps);

	// Fill in the values
	for (unsigned long i = 0; i < NewSteps+1; i++)
	{
		for (unsigned long j = 0; j < Steps+1; j++)
		{
			double Integral = SubIntegral(sLow+i*sStep,vLow+j*vStep,NewTimeStep,Proc,Old);
			Result.SetValue(i,j,Integral);
		}
	}

	return Result;
}

#endif
