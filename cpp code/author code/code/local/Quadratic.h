#ifndef QUADRATIC_H
#define QUADRATIC_H

#include "LVProcess.h" 
#include "MersenneTwister.h"
#include "MiscMaths.h"
#include <math.h>

class	QuadraticProcess : public LVProcess
{
public:
	QuadraticProcess(double Spot_, double Rate_, double Vol_, double b_);
	void Evolve(double Time, MTRand& MT); 
	void Evolve(double Time, unsigned long N, MTRand& MT);
	void EvolveMilstein(double Time, MTRand& MT); 
	void EvolveMilstein(double Time, unsigned long N, MTRand& MT); 
	double InitVol() const;
	double InitVol(double Base) const;
	double InitLow(double Time, double Deviations) const;
	double InitHigh(double Time, double Deviations) const;
	double InitLow(double Base, double Time, double Deviations) const;
	double InitHigh(double Base, double Time, double Deviations) const;
	double Density(double Time, double Value) const;
	double Density(double Time, double Value, double Base) const;
	double ApprDensity(double Time, double Value, double Base) const;
	double AS1Density(double Time, double Value, double Base) const;
	double GetValue() const;
	double GetRate() const;
	void SetValue(double NewValue); 
private:
	double Spot;
	double Rate;
	double Vol;
	double b;
	double Volb,K1, K2, K3;
};


double QuadraticProcess::GetValue() const
{
	return Spot;
}

double QuadraticProcess::GetRate() const
{
	return Rate;
}

QuadraticProcess::QuadraticProcess(double Spot_, double Rate_, double Vol_, double b_) : 
Spot(Spot_), Rate(Rate_), Vol(Vol_), b(b_)
{
	Volb = Vol*b;
	K1 = 4*Rate*Rate*Volb - 4*Rate*Volb*Volb*Volb + Volb*Volb*Volb*Volb*Volb;
	K2 = -8*Rate*Vol*(Rate+Volb*Volb);
	K3 = 2*Rate*Rate*Vol*Vol;
}

void QuadraticProcess::Evolve(double Time, MTRand& MT)
{
	Spot += Rate*Spot*Time + Vol*Spot*(Spot+b)*sqrt(Time)*MT.randNorm();
}

void QuadraticProcess::EvolveMilstein(double Time, MTRand& MT)
{
	double R = sqrt(Time)*MT.randNorm();
	Spot += Rate*Spot*Time + Vol*Spot*(Spot+b)*R + 0.5*Vol*Vol*Spot*(2*Spot+b)*(Spot+b)*(R*R-Time);
}


void QuadraticProcess::Evolve(double Time, unsigned long N, MTRand& MT)
{
	double TimeStep = Time/static_cast<double>(N); 
	unsigned long i = 0; 

	while(i < N && Spot > 0)
	{
		Evolve(TimeStep,MT);
		i++; 
	}
}

void QuadraticProcess::EvolveMilstein(double Time, unsigned long N, MTRand& MT)
{
	double TimeStep = Time/static_cast<double>(N); 
	unsigned long i = 0; 

	while(i < N && Spot > 0)
	{
		EvolveMilstein(TimeStep,MT);
		i++; 
	}
}

double QuadraticProcess::InitVol() const
{
	return Vol*(Spot+b);
}

double QuadraticProcess::InitVol(double Base) const
{
	return Vol*(Base+b);
}

double QuadraticProcess::InitLow(double Time, double Deviations) const
{
	return InitLow(Spot,Time,Deviations);
}

double QuadraticProcess::InitHigh(double Time, double Deviations) const
{
	return InitHigh(Spot,Time,Deviations);
}

double QuadraticProcess::InitLow(double Base, double Time, double Deviations) const
{
	double IV = InitVol(Base);
	return Base*exp((Rate-0.5*IV*IV)*Time - Deviations*IV*sqrt(Time));
}

double QuadraticProcess::InitHigh(double Base, double Time, double Deviations) const
{
	double IV = InitVol(Base);
	return Base*exp((Rate-0.5*IV*IV)*Time + Deviations*IV*sqrt(Time));
}

double QuadraticProcess::Density(double Time, double Value) const
{
	return Density(Time,Value,Spot);
}

double QuadraticProcess::Density(double Time, double Value, double Base) const
{
	return AS1Density(Time,Value,Base);
}


void QuadraticProcess::SetValue(double NewValue) 
{
	Spot = NewValue;
}

double QuadraticProcess::ApprDensity(double Time, double Value, double Base) const
{
	if ( Value == Base )
	{
		// Need to compute the value as a limit in this case
		return (1 - 0.125*pow(Volb,-3.0)*(K1 + Volb*exp(Volb*Value)*K2 + 2*Volb*exp(2.0*Volb*Value)*K3)*Time)/(Vol*Value*(Value+b)*sqrt(2*PI*Time));
	}
	else
	{
		double y = log(Value/(Vol*Value+Volb))/Volb;
		double y0 = log(Base/(Vol*Base+Volb))/Volb;
		double E = exp(Volb*y);
		double E0 = exp(Volb*y0);
		double Result = 1/(Vol*Value*(Value+b)*sqrt(2*PI*Time));
		Result *= (1 - 0.125*pow(Volb,-3.0)*(K1 + K2*(E-E0)/(Value-Base) + K3*(E*E-E0*E0)/(Value-Base))*Time);
		Result *= exp(-(y-y0)*(y-y0)/2/Time + (Rate/Volb - Volb/2)*(y-y0) - Rate*Vol/Volb/Volb*(E-E0) + log((Vol*E-1)/(Vol*E0-1)));
		return Result;
	}
}

double QuadraticProcess::AS1Density(double Time, double Value, double Base) const
{
	if (Value == Base)
	{
		return AS1Density(Time,Value,Value+0.0001);
	}
	else
	{
		double c = Vol*b;
		double y = 1/c*log(Value/(Value+b));
		double yInit = 1/c*log(Base/(Base+b));
		double yDiff = y - yInit;
		double eDiff = exp(y*c) - exp(yInit*c);
		double e2Diff = exp(2*y*c) - exp(2*yInit*c);

		double Factor1 = 1 - Time*((4*Rate*Rate*c + pow(c,5) - 4*Rate*c*c*c)*yDiff - 8*Rate*(Rate+c*c)*eDiff + 2*Rate*Rate*e2Diff)/8.0/yDiff/c/c/c;
		double Factor2 = (exp(y*c)-1)/(exp(yInit*c)-1);
		double Factor3 = exp(-yDiff*yDiff/2/Time + (Rate/c - c/2.0)*yDiff - Rate/c/c*eDiff);
		return Factor1*Factor2*Factor3/sqrt(2*PI*Time)/Vol/Value/(Value+b);
	}
}


#endif
