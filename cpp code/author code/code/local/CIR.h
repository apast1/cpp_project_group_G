#ifndef CIR_H
#define CIR_H

#include "MersenneTwister.h"
#include "MiscMaths.h"
#include <math.h>
#include <boost/math/special_functions/bessel.hpp>

using boost::math::cyl_bessel_i;

class CIRProcess 
{
public:
	CIRProcess(double Rate_, double theta_, double kappa_, double Vol_);
	void Evolve(double Time, MTRand& MT); 
	void Evolve(double Time, unsigned long N, MTRand& MT);
	void EvolveMilstein(double Time, MTRand& MT); 
	void EvolveMilstein(double Time, unsigned long N, MTRand& MT); 
	void CleverEvolve(double Time, unsigned long MinSteps, double Limit, MTRand& MT);
	double InitVol() const;
	double InitLow(double Time, double Deviations) const;
	double InitHigh(double Time, double Deviations) const;
	double Density(double Time, double Value) const;
	double Density(double Time, double Value, double Base) const;
	double ApprZeroDensity(double Time, double Value, double Base) const;
	double ApprZeroDensity(double Time, double Value) const;
	double ExactDensity(double Time, double Value, double Base) const;
	double ExactDensity(double Time, double Value) const;
	double ZeroProb(double Time, unsigned long N) const;
	double GetValue() const;
	void SetValue(double NewValue); 
private:
	double Rate;
	double theta;
	double kappa;
	double Vol;
};


double CIRProcess::GetValue() const
{
	return Rate;
}


CIRProcess::CIRProcess(double Rate_, double theta_, double kappa_, double Vol_) : 
Rate(Rate_), theta(theta_), kappa(kappa_), Vol(Vol_)
{
}

void CIRProcess::Evolve(double Time, MTRand& MT)
{
	if (Rate != 0)
	{
		Rate += kappa*(theta-Rate)*Time + Vol*sqrt(Rate*Time)*MT.randNorm(); 
		if (Rate <= 0)
		{
			Rate = 0;
		}
	}	
}

void CIRProcess::EvolveMilstein(double Time, MTRand& MT)
{
	if (Rate != 0)
	{
		double z = MT.randNorm(); 
		Rate += kappa*(theta-Rate)*Time + Vol*sqrt(Rate*Time)*z + 0.25*Vol*Vol*(z*z-1)*Time; 
		if (Rate <= 0)
		{
			Rate = 0;
		}
	}	
}

/*
void CIRProcess::CleverEvolve(double Time, unsigned long MinSteps, double Limit, MTRand& MT)
{
	double TimeLapsed = 0; 
	while (TimeLapsed < Time)
	{
		double IVar = pow(Vol*pow(Rate,beta-1),2.0); 
		double AltTime = Limit*Limit/IVar;
		double TestTimeStep = min(Time/static_cast<double>(MinSteps),AltTime); 
		if (TestTimeStep < 0.0001)
		{
			Rate = 0;
			TimeLapsed = Time; 
		}
		else
		{
			double TimeStep = min(Time-TimeLapsed,TestTimeStep); 
			Evolve(TimeStep,MT);
			TimeLapsed += TimeStep;
		}
	}
}
*/

void CIRProcess::Evolve(double Time, unsigned long N, MTRand& MT)
{
	double TimeStep = Time/static_cast<double>(N); 
	unsigned long i = 0; 

	while(i < N && Rate > 0)
	{
		Evolve(TimeStep,MT);
		i++; 
	}
}

void CIRProcess::EvolveMilstein(double Time, unsigned long N, MTRand& MT)
{
	double TimeStep = Time/static_cast<double>(N); 
	unsigned long i = 0; 

	while(i < N && Rate > 0)
	{
		EvolveMilstein(TimeStep,MT);
		i++; 
	}
}

double CIRProcess::InitVol() const
{
	return Vol/sqrt(Rate); 
}

double CIRProcess::InitLow(double Time, double Deviations) const
{
	double IV = InitVol(); 
	return Rate*exp((Rate-0.5*IV*IV)*Time - IV*sqrt(Time)*Deviations); 
}

double CIRProcess::InitHigh(double Time, double Deviations) const
{
	double IV = InitVol();
	return Rate*exp((Rate-0.5*IV*IV)*Time + IV*sqrt(Time)*Deviations); 
}

double CIRProcess::ExactDensity(double Time, double Value) const
{
	return ExactDensity(Time,Value,Rate); 
}

double CIRProcess::ExactDensity(double Time, double Value, double Base) const
{
	double c = 2*kappa/(Vol*Vol*(1-exp(-kappa*Time))); 
	double q = 2*kappa*theta/(Vol*Vol)-1; 
	double alpha = c*Base*exp(-kappa*Time); 
	double nu = c*Rate;
	return c*exp(-alpha-nu+kappa*Time*q/2)*cyl_bessel_i(q,2*sqrt(alpha*nu)); 
}

double CIRProcess::Density(double Time, double Value) const
{
	return Density(Time,Value,Rate);
}

double CIRProcess::Density(double Time, double Value, double Base) const
{
	double root = sqrt(2.0); 
	double APrime = root*kappa*theta/Vol-Vol/(2*root); 
	double BPrime = root*kappa/Vol;
	double A = 2.0*root*APrime/Vol;
	double B = BPrime*Vol/(2.0*root);
	double Var = Vol*Vol; 
	double Result = 1/(Vol*sqrt(2.0*PI*Time*Value))*exp(-0.25*A*log(Base/Value)+2.0*B*(Base-Value)/Var);
	Result *= exp(-2*(Base-Value)*(Base-Value)/(Time*Var)); 
	double Mess = (APrime/sqrt(Value)-BPrime*sqrt(Value));
	Result *= (1 - 0.25*Mess*Time + 0.5*Mess*Mess*Time*Time); 
	return Result; 
}

void CIRProcess::SetValue(double NewValue) 
{
	Rate = NewValue;
}

#endif