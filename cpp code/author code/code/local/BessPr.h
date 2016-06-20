#ifndef BESSEL_PROCESS_H
#define BESSEL_PROCESS_H

#include "MersenneTwister.h"

class BesselProcess
{
public:
	BesselProcess(double Spot_, double A_);
	void Evolve(double TimeStep, MTRand& MT);
	void Evolve(double Time, unsigned long N, MTRand& MT);
	double GetValue() const;
	void SetValue(double NewValue);
private:
	double Spot;
	double A;
};

BesselProcess::BesselProcess(double Spot_, double A_) 
{
	Spot = Spot_;
	A = A_;
}

void BesselProcess::Evolve(double TimeStep, MTRand& MT)
{
	if (Spot != 0)
	{
		Spot += (A/Spot)*TimeStep + sqrt(TimeStep)*MT(); 
	}
}

void BesselProcess::Evolve(double Time, unsigned long N, MTRand& MT)
{
	double TimeStep = Time/static_cast<double>(N); 
	for (unsigned long i = 0; i < N; i++)
	{
		Evolve(TimeStep,MT);
	}
}

double BesselProcess::GetValue() const
{
	return Spot;
}

void BesselProcess::SetValue(double NewValue)
{
	Spot = NewValue;
}


#endif