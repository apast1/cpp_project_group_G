#ifndef BSPROCESS_H
#define BSPROCESS_H


#include <complex> 
#include <cmath>
#include <vector>

#include "MiscMaths.h"
#include "Constants.h"
#include "MersenneTwister.h"
#include "LVProcess.h" 

using std::vector; 

double BSLogDensity(double logInit, double logS, double Rate, double Vol, double Time);


class BSProcess : public LVProcess
{
public: 
	BSProcess(double Spot_, double Rate_, double Vol_); 
	double GetVol() const;
	double GetRate() const;
	double GetValue() const;
	double Density(double Time, double Value) const;
	double Density(double Time, double Value, double Base) const;
	double EuroCall(double Strike, double Time) const; 
	double EuroPut(double Strike, double Time) const; 
	double InitLow(double Time, double Deviations) const;
	double InitHigh(double Time, double Deviations) const;
	double InitLow(double Base, double Time, double Deviations) const;
	double InitHigh(double Base, double Time, double Deviations) const;
	void Evolve(double Time, MTRand& MT); 
	void Evolve(double Time, unsigned long N, MTRand& MT);
	void SetValue(double NewValue); 
	Complex Char(double Time, Complex z) const;
private:
	double Spot; 
	double Rate; 
	double Vol;
};

double BSProcess::InitLow(double Base, double Time, double Deviations) const
{
	return Base*exp((Rate-0.5*Vol*Vol)*Time - Deviations*Vol*sqrt(Time)); 
}

double BSProcess::InitHigh(double Base, double Time, double Deviations) const
{
	return Base*exp((Rate-0.5*Vol*Vol)*Time + Deviations*Vol*sqrt(Time)); 
}

double BSProcess::InitLow(double Time, double Deviations) const
{
	return Spot*exp((Rate-0.5*Vol*Vol)*Time - Deviations*Vol*sqrt(Time)); 
}

double BSProcess::InitHigh(double Time, double Deviations) const
{
	return Spot*exp((Rate-0.5*Vol*Vol)*Time + Deviations*Vol*sqrt(Time)); 
}

BSProcess::BSProcess(double Spot_, double Rate_, double Vol_) : Spot(Spot_), Rate(Rate_), Vol(Vol_)
{
}

void BSProcess::Evolve(double Time, MTRand& MT)
{
	Spot *= exp((Rate-0.5*Vol*Vol)*Time + Vol*sqrt(Time)*MT()); 
}

void BSProcess::Evolve(double Time, unsigned long N, MTRand& MT)
{
	double TimeStep = Time/static_cast<double>(N); 
	for (unsigned long i = 0; i < 0; i++)
	{
		Spot += Rate*Spot*TimeStep + Vol*Spot*sqrt(TimeStep)*MT(); 
	}
}


double BSProcess::Density(double Time, double Value) const
{
	return Density(Time,Value,Spot); 
}

double BSProcess::Density(double Time, double Value, double Base) const
{
	return exp(-pow(log(Value/Base)-(Rate-0.5*Vol*Vol)*Time,2.0)/(2.0*Vol*Vol*Time))/(Vol*Value*sqrt(2.0*PI*Time));
//	return 1;

}

double BSProcess::GetValue() const
{
	return Spot;
}

double BSProcess::GetVol() const
{
	return Vol;
}


Complex BSProcess::Char(double Time, Complex z) const
{
	// The characteristic function for Black-Scholes process
	Complex y; 
	Complex I(0.0,1.0);
	y = (log(Spot)+ (Rate-0.5*Vol*Vol)*Time)*I*z-0.5*Vol*Vol*Time*z*z;
	return exp(y); 
}

double BSProcess::EuroCall(double Strike, double Time) const
{
	// Returns the price of European call using closed-form formulae
	double d = (log(Spot/Strike)+(Rate+0.5*Vol*Vol)*Time)/(Vol*sqrt(Time)); 
	double e = d-Vol*sqrt(Time); 

	return (Spot*Normal(d) - exp(-Rate*Time)*Strike*Normal(e)); 
}

double BSProcess::EuroPut(double Strike, double Time) const
{
	// Returns the price of European put using closed-form formulae
	double d = (log(Spot/Strike)+(Rate+0.5*Vol*Vol)*Time)/(Vol*sqrt(Time));
	double e = d - Vol*sqrt(Time);

	return Strike*exp(-Rate*Time)*Normal(-e) - Spot*Normal(-d); 
}


double BSProcess::GetRate() const
{
	return Rate;
}

void BSProcess::SetValue(double NewValue)
{
	Spot = NewValue;
}


#endif
