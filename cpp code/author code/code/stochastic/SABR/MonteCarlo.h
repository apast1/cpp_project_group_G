#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <vector>

#include "MersenneTwister.h"
#include "LVProcess.h"

using std::vector;

double MCLookbackCall(double Strike, vector<double> Times, LVProcess& Proc, MTRand& MT, unsigned long N);

vector<double> MCVector(vector<double> Times, unsigned long N, LVProcess& Proc, MTRand& MT);

double MaxEl(vector<double> Vec);

double MCLookbackCall(double Strike, vector<double> Times, LVProcess& Proc, MTRand& MT, unsigned long N)
{
	double Sum = 0;
	double InitSpot = Proc.GetValue(); 

	for (unsigned long i = 0; i < N*N; i++)
	{
		vector<double> Vec = MCVector(Times,N,Proc,MT); 
		Sum += max(MaxEl(Vec)-Strike,0.0);
	}

	return exp(-Proc.GetRate()*Times.back())*Sum/static_cast<double>(N*N); 
}


vector<double> MCVector(vector<double> Times, unsigned long N, LVProcess& Proc, MTRand& MT)
{
	vector<double> Result;
	double InitSpot = Proc.GetValue(); 
	for (unsigned long i = 0; i < Times.size(); i++)
	{
		double TimeStep = 0;
		if (i == 0)
		{
			TimeStep = Times[0];
		}
		else
		{
			TimeStep = Times[i] - Times[i-1]; 
		}
		Proc.Evolve(TimeStep,N,MT);
		Result.push_back(Proc.GetValue()); 
	}
	Proc.SetValue(InitSpot); 

	return Result;
}

double MaxEl(vector<double> Vec)
{
	double Result = Vec[0];
	for (unsigned long i = 1; i < Vec.size(); i++)
	{
		if (Vec[i] > Result)
		{
			Result = Vec[i];
		}
	}
	return Result;
}

#endif