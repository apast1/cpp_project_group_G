#ifndef CEV_H
#define CEV_H

#include "LVProcess.h" 
#include "MersenneTwister.h"
#include "MiscMaths.h"
#include <math.h>


class CEVProcess : public LVProcess
{
public:
	CEVProcess(double Spot_, double beta_, double Rate_, double Vol_);
	void Evolve(double Time, MTRand& MT); 
	void Evolve(double Time, unsigned long N, MTRand& MT);
	void EvolveMilstein(double Time, MTRand& MT); 
	void EvolveMilstein(double Time, unsigned long N, MTRand& MT); 
	void CleverEvolve(double Time, unsigned long MinSteps, double Limit, MTRand& MT);
	double InitVol() const;
	double InitVol(double Base) const;
	double EuroPut(double Time, double Strike, unsigned long Steps) const; 
	double InitLow(double Time, double Deviations) const;
	double InitHigh(double Time, double Deviations) const;
	double InitLow(double Base, double Time, double Deviations) const;
	double InitHigh(double Base, double Time, double Deviations) const;
	double EuroPutViaAppr(double Time, double Strike, unsigned long Steps) const; 
	double EuroCallViaAppr(double Time, double Strike, unsigned long Steps) const;
	double Density(double Time, double Value) const;
	double Density(double Time, double Value, double Base) const;
	double HLDensity(double Time, double Value, double Base) const;
	double ApprZeroDensity(double Time, double Value, double Base) const;
	double ApprZeroDensity(double Time, double Value) const;
	double ExactDensity(double Time, double Value, double Base) const;
	double ExactDensity(double Time, double Value) const;
	double ZeroProb(double Time, unsigned long N) const;
	double ASDensity(double Time, double Value, double Base) const;
	double AS1Density(double Time, double Value, double Base) const;
	double GetValue() const;
	double GetRate() const;
	void SetValue(double NewValue); 
private:
	double Spot;
	double beta;
	double Rate;
	double Vol;
};


double CEVProcess::GetValue() const
{
	return Spot;
}

double CEVProcess::GetRate() const
{
	return Rate;
}

CEVProcess::CEVProcess(double Spot_, double beta_, double Rate_, double Vol_) : 
Spot(Spot_), beta(beta_), Rate(Rate_), Vol(Vol_)
{
}

void CEVProcess::Evolve(double Time, MTRand& MT)
{
	if (Spot != 0)
	{
		Spot += Rate*Spot*Time + Vol*pow(Spot,beta)*sqrt(Time)*MT.randNorm(); 
		if (Spot <= 0)
		{
			Spot = 0;
		}
	}	
}

void CEVProcess::EvolveMilstein(double Time, MTRand& MT)
{
	if (Spot != 0)
	{
		double z = MT.randNorm(); 
		Spot += Rate*Spot*Time + Vol*pow(Spot,beta)*sqrt(Time)*z + 0.5*beta*Vol*Vol*pow(Spot,2*beta-1)*Time*(z*z-1); 
		if (Spot <= 0)
		{
			Spot = 0;
		}
	}	
}

void CEVProcess::CleverEvolve(double Time, unsigned long MinSteps, double Limit, MTRand& MT)
{
	double TimeLapsed = 0; 
	while (TimeLapsed < Time)
	{
		double IVar = pow(Vol*pow(Spot,beta-1),2.0); 
		double AltTime = Limit*Limit/IVar;
		double TestTimeStep = min(Time/static_cast<double>(MinSteps),AltTime); 
		if (TestTimeStep < 0.0001)
		{
			Spot = 0;
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

void CEVProcess::Evolve(double Time, unsigned long N, MTRand& MT)
{
	double TimeStep = Time/static_cast<double>(N); 
	unsigned long i = 0; 

	while(i < N && Spot > 0)
	{
		Evolve(TimeStep,MT);
		i++; 
	}
}

void CEVProcess::EvolveMilstein(double Time, unsigned long N, MTRand& MT)
{
	double TimeStep = Time/static_cast<double>(N); 
	unsigned long i = 0; 

	while(i < N && Spot > 0)
	{
		EvolveMilstein(TimeStep,MT);
		i++; 
	}
}

double CEVProcess::InitVol() const
{
	return Vol*pow(Spot,beta-1); 
}

double CEVProcess::InitVol(double Base) const
{
	return Vol*pow(Base,beta-1); 
}

double CEVProcess::InitLow(double Time, double Deviations) const
{
	double IV = InitVol(); 
	return Spot*exp((Rate-0.5*IV*IV)*Time - IV*sqrt(Time)*Deviations); 
}

double CEVProcess::InitHigh(double Time, double Deviations) const
{
	double IV = InitVol();
	return Spot*exp((Rate-0.5*IV*IV)*Time + IV*sqrt(Time)*Deviations); 
}

double CEVProcess::InitLow(double Base, double Time, double Deviations) const
{
	double IV = InitVol(Base); 
	return Base*exp((Rate-0.5*IV*IV)*Time - IV*sqrt(Time)*Deviations); 
}

double CEVProcess::InitHigh(double Base, double Time, double Deviations) const
{
	double IV = InitVol(Base);
	return Base*exp((Rate-0.5*IV*IV)*Time + IV*sqrt(Time)*Deviations); 
}
 

double CEVProcess::ApprZeroDensity(double Time, double Value) const
{
	return ApprZeroDensity(Time,Value,Spot);
}

double CEVProcess::ApprZeroDensity(double Time, double Value, double Base) const
{
	double LocalTime = Vol*Vol*Time;
	double fav = (Value + Base)/2.0; 
	double Result = 1.0; 
	Result += 0.125*pow(fav,2*(beta-1.0))*beta*(beta-2.0)*LocalTime; 
	Result += beta*(beta-2.0)*(3*beta-2.0)*(3.0*beta-4.0)*pow(fav,4.0*beta-4.0)*LocalTime*LocalTime/128.0;
	Result *= pow(Value,-beta)*pow(Base/Value,beta/2.0)*exp(-(pow(pow(Value,1-beta)-pow(Base,1-beta),2.0))/(2*LocalTime*(1-beta)*(1-beta)))/(sqrt(2*PI*LocalTime)); 
	return Result;
}

double CEVProcess::Density(double Time, double Value) const
{
	return Density(Time,Value,Spot);
}

double CEVProcess::Density(double Time, double Value, double Base) const
{
	return AS1Density(Time,Value,Base);
/*	double NewTime = Time;
	if ( Rate != 0)
	{
		NewTime = (exp(2*Rate*Time*(beta-1)) - 1)/(2.0*Rate*(beta-1)); 
	}
	return exp(-Rate*Time)*ApprZeroDensity(NewTime,exp(-Rate*Time)*Value,Base); 
*/
}

double CEVProcess::HLDensity(double Time, double Value, double Base) const
{
	double NewTime = Time;
	if ( Rate != 0)
	{
		NewTime = (exp(2*Rate*Time*(beta-1)) - 1)/(2.0*Rate*(beta-1)); 
	}
	return exp(-Rate*Time)*ApprZeroDensity(NewTime,exp(-Rate*Time)*Value,Base); 
}

double CEVProcess::EuroPut(double Time, double Strike, unsigned long Steps) const
{
	double IV = InitVol(); 
	double Term = exp((Rate-0.5*IV*IV)*Time); 
	double LowerLim = Spot*Term*exp(-50.0*IV*sqrt(Time)); 
	double UpperLim = min(Strike,Spot*Term*exp(40.0*IV*sqrt(Time)));
	double Step = (UpperLim-LowerLim)/static_cast<double>(Steps); 
	double preResult = max(0.0,Strike-LowerLim)*Density(Time,LowerLim);
	preResult += max(0.0,Strike-UpperLim)*Density(Time,UpperLim); 
	double pneq0 = Density(Time,LowerLim)+Density(Time,UpperLim);

	for (unsigned long i = 1; i < 2*(Steps-1); i++)
	{
		double Value = LowerLim+i*Step; 
		double PayOff = max(0.0,Strike-Value); 
		double DRN = Density(Time,Value); 
		double Coeff;
		if (i > 0 && i%2)
		{
			Coeff = 4.0; 
		}
		else
		{
			Coeff = 2.0;
		}
		pneq0 += Coeff*DRN;
		preResult += Coeff*DRN*PayOff;
	}	
	preResult *= Step/3.0;
	pneq0 *= Step/3.0;
	return (1-pneq0)*Strike+preResult; 
}

double CEVProcess::EuroPutViaAppr(double Time, double Strike, unsigned long Steps) const
{
	double IV = InitVol(); 
	double Term = exp((Rate-0.5*IV*IV)*Time); 
	double LowerLim = Spot*Term*exp(-50.0*IV*sqrt(Time)); 
	double UpperLim = min(Strike,Spot*Term*exp(40.0*IV*sqrt(Time)));
	double Step = (UpperLim-LowerLim)/static_cast<double>(Steps); 
	double preResult = max(0.0,Strike-LowerLim)*Density(Time,LowerLim);
	preResult += max(0.0,Strike-UpperLim)*Density(Time,UpperLim); 
	double pneq0 = Density(Time,LowerLim)+ Density(Time,UpperLim);

	for (unsigned long i = 1; i < 2*(Steps-1); i++)
	{
		double Value = LowerLim+i*Step; 
		double PayOff = max(0.0,Strike-Value); 
		double DRN = Density(Time,Value); 
		double Coeff;
		if (i > 0 && i%2)
		{
			Coeff = 4.0; 
		}
		else
		{
			Coeff = 2.0;
		}
		pneq0 += Coeff*DRN;
		preResult += Coeff*DRN*PayOff;
	}	
	preResult *= Step/3.0;
	pneq0 *= Step/3.0;
	return (1-pneq0)*Strike+preResult; 
}

double CEVProcess::EuroCallViaAppr(double Time, double Strike, unsigned long Steps) const
{
	return Spot - exp(-Rate*Time)*Strike + EuroPutViaAppr(Time,Strike,Steps); 
}

double CEVProcess::ZeroProb(double Time, unsigned long Steps) const
{
	double IV = InitVol(); 
	double Term = exp((Rate-0.5*IV*IV)*Time); 
	double LowerLim = Spot*Term*exp(-50.0*IV*sqrt(Time)); 
	double UpperLim = Spot*Term*exp(10.0*IV*sqrt(Time));
	double Step = (UpperLim-LowerLim)/static_cast<double>(2*Steps+1);

	double pneq0 = Density(Time,LowerLim)+Density(Time,UpperLim);

	for (unsigned long i = 1; i < 2*(Steps-1); i++)
	{
		double Value = LowerLim+i*Step; 
		double DRN = Density(Time,Value); 
		double Coeff;
		if (i%2)
		{
			Coeff = 4.0; 
		}
		else
		{
			Coeff = 2.0;
		}
		pneq0 += Coeff*DRN;
	}	
	pneq0 *= Step/3.0;
	return 1.0-pneq0; 
}

void CEVProcess::SetValue(double NewValue) 
{
	Spot = NewValue;
}

double CEVProcess::AS1Density(double Time, double Value, double Base) const
{
	if (Value == Base)
	{
		return AS1Density(Time,Value,Base+0.0001);
	}
	else
	{
		double r = Rate;
	 	double b = beta; 
		double y = pow(Value,1-beta)/Vol/(1-beta);
		double yInit = pow(Base,1-beta)/Vol/(1-beta);
		return (1.0+1/(y-yInit)*(-(24.0*r*r*y*y*y*y*b*b+4.0*r*r*y*y*y*y-16.0*r*r*y*y*y*y*b+4.0*r*r*y*y*y*y*b*b*b*b-16.0*r*r*y*y*y*y*b*b*b+60.0*r*y*y*b*b-48.0*r*y*y*b-24.0*r*b*b*b*y*y+3.0*b*b+12.0*r*y*y-6.0*b)/pow(-1.0+b,2.0)/y/24.0+(24.0*r*r*yInit*yInit*yInit*yInit*b*b+4.0*r*r*yInit*yInit*yInit*yInit-16.0*r*r*yInit*yInit*yInit*yInit*b+4.0*r*r*yInit*yInit*yInit*yInit*b*b*b*b-16.0*r*r*yInit*yInit*yInit*yInit*b*b*b+60.0*r*yInit*yInit*b*b-48.0*r*yInit*yInit*b-24.0*r*b*b*b*yInit*yInit+3.0*b*b+12.0*r*yInit*yInit-6.0*b)/pow(-1.0+b,2.0)/yInit/24.0)*Time)*sqrt(2.0)/sqrt(3.141592653589793*Time)*exp(-pow(y-yInit,2.0)/Time/2.0)*exp(-1/(-1.0+b)*r*y*y/2.0+1/(-1.0+b)*r*y*y*b-1/(-1.0+b)*r*y*y*b*b/2.0+1/(-1.0+b)*b*log(y)/2.0+1/(-1.0+b)*r*yInit*yInit/2.0-1/(-1.0+b)*r*yInit*yInit*b+1/(-1.0+b)*r*yInit*yInit*b*b/2.0-1/(-1.0+b)*b*log(yInit)/2.0)/2.0/Vol/pow(Value,beta);
	}
}

double CEVProcess::ASDensity(double Time, double Value, double Base) const
{
	if (Value == Base)
	{
		return ASDensity(Time,Value,Base+0.0001);
	}
	else
	{
		double r = Rate;
	 	double b = beta; 
		double y = pow(Value,1-beta)/Vol/(1-beta);
		double yInit = pow(Base,1-beta)/Vol/(1-beta);
		
		double MapleGenVar2 = 1.0/2.0;
		double MapleGenVar5 = 1.0;      
		double MapleGenVar7 = 1/(y-yInit)*(-(24.0*r*r*y*y*y*y*b*b+4.0*r*r*y*y*y*y-16.0*r*r*y*y*y*y*b+4.0*r*r*y*y*y*y*b*b*b*b-16.0*r*r*y*y*y*y*b*b*b+60.0*r*y*y*b*b-48.0*r*y*y*b-24.0*r*b*b*b*y*y+3.0*b*b+12.0*r*y*y-6.0*b)/pow(-1.0+b,2.0)/y/24.0+(24.0*r*r*yInit*yInit*yInit*yInit*b*b+4.0*r*r*yInit*yInit*yInit*yInit-16.0*r*r*yInit*yInit*yInit*yInit*b+4.0*r*r*yInit*yInit*yInit*yInit*b*b*b*b-16.0*r*r*yInit*yInit*yInit*yInit*b*b*b+60.0*r*yInit*yInit*b*b-48.0*r*yInit*yInit*b-24.0*r*b*b*b*yInit*yInit+3.0*b*b+12.0*r*yInit*yInit-6.0*b)/pow(-1.0+b,2.0)/yInit/24.0)*Time;      
		double MapleGenVar9 = 1/(pow(y-yInit,2.0));      
		double MapleGenVar15 = 11.0/32.0*yInit*b*b-yInit*b/8.0-r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit/36.0-r*r*r*yInit*yInit*yInit*yInit*y*y*y/12.0+7.0/8.0*r*b*b*b*y*y*y+r*r*y*y*y*y*y*b*b*b*b*b/8.0-9.0/16.0*b*b*b*b*y*y*y*r-3.0/16.0*r*r*y*y*y*y*y*b*b+r*r*y*y*y*y*y*b*b*b/3.0-yInit*yInit*r*r*y*y*y/12.0-r*r*r*y*y*y*y*y*yInit*yInit/12.0-7.0/24.0*r*r*y*y*y*y*y*b*b*b*b-9.0/16.0*r*y*y*y*b*b+r*r*y*y*y*y*y*b/24.0+b*b*b*b*b*y*y*y*r/8.0-r*r*y*y*y*y*y*b*b*b*b*b*b/48.0+r*y*y*y*b/8.0+b*y/4.0-9.0/32.0*yInit*b*b*b+9.0/128.0*yInit*b*b*b*b+23.0/3.0*yInit*yInit*r*r*b*b*b*y*y*y;      
		double MapleGenVar14 = MapleGenVar15+4.0*r*r*b*b*b*b*b*yInit*yInit*y*y*y-13.0/12.0*r*r*r*b*b*b*b*b*b*yInit*yInit*yInit*yInit*y*y*y-4.0*yInit*yInit*r*r*b*b*y*y*y+3.0*r*r*r*b*b*b*b*b*yInit*yInit*yInit*yInit*y*y*y-9.0/4.0*r*r*r*yInit*yInit*yInit*yInit*b*b*y*y*y+yInit*yInit*r*r*b*y*y*y-31.0/4.0*yInit*yInit*r*r*b*b*b*b*y*y*y+r*r*r*b*b*b*b*b*b*b*yInit*yInit*yInit*yInit*y*y*y/6.0+2.0/3.0*r*r*r*yInit*yInit*yInit*yInit*b*y*y*y+25.0/6.0*r*r*r*b*b*b*yInit*yInit*yInit*yInit*y*y*y-55.0/12.0*r*r*r*b*b*b*b*yInit*yInit*yInit*yInit*y*y*y-5.0/6.0*r*r*b*b*b*b*b*b*yInit*yInit*y*y*y-13.0/12.0*r*r*r*y*y*y*y*y*b*b*b*b*b*b*yInit*yInit+2.0/9.0*r*r*r*r*y*y*y*y*y*b*b*b*b*b*b*b*yInit*yInit*yInit*yInit+2.0/9.0*r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit*b-35.0/18.0*r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit*b*b*b*b+14.0/9.0*r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit*b*b*b-9.0/4.0*r*r*r*y*y*y*y*y*yInit*yInit*b*b+2.0/3.0*r*r*r*y*y*y*y*y*yInit*yInit*b+25.0/6.0*r*r*r*y*y*y*y*y*b*b*b*yInit*yInit-7.0/9.0*r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit*b*b-55.0/12.0*r*r*r*y*y*y*y*y*b*b*b*b*yInit*yInit;      
		MapleGenVar15 = MapleGenVar14+14.0/9.0*r*r*r*r*y*y*y*y*y*b*b*b*b*b*yInit*yInit*yInit*yInit+3.0*r*r*r*y*y*y*y*y*b*b*b*b*b*yInit*yInit-7.0/9.0*r*r*r*r*y*y*y*y*y*b*b*b*b*b*b*yInit*yInit*yInit*yInit-r*r*r*r*y*y*y*y*y*b*b*b*b*b*b*b*b*yInit*yInit*yInit*yInit/36.0+yInit*r*r*y*y*y*y/24.0+35.0/16.0*yInit*r*r*y*y*y*y*b*b+25.0/6.0*yInit*r*r*y*y*y*y*b*b*b*b-25.0/6.0*yInit*r*r*y*y*y*y*b*b*b-17.0/8.0*yInit*r*r*y*y*y*y*b*b*b*b*b+7.0/16.0*yInit*r*r*y*y*y*y*b*b*b*b*b*b-13.0/24.0*yInit*r*r*y*y*y*y*b+r*r*r*r*yInit*y*y*y*y*y*y*y*y/72.0-7.0/9.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b*b+r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b*b*b*b*b/72.0-r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b*b*b*b/9.0+35.0/36.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b-7.0/9.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b+7.0/18.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b+7.0/18.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b*b*b-r*r*r*r*yInit*y*y*y*y*y*y*y*y*b/9.0+r*r*r*yInit*y*y*y*y*y*y/12.0;      
		double MapleGenVar13 = MapleGenVar15+13.0/12.0*r*r*r*yInit*y*y*y*y*y*y*b*b*b*b*b*b+9.0/4.0*r*r*r*yInit*y*y*y*y*y*y*b*b-25.0/6.0*r*r*r*yInit*y*y*y*y*y*y*b*b*b+55.0/12.0*r*r*r*yInit*y*y*y*y*y*y*b*b*b*b-3.0*r*r*r*yInit*y*y*y*y*y*y*b*b*b*b*b-r*r*r*yInit*y*y*y*y*y*y*b*b*b*b*b*b*b/6.0-2.0/3.0*r*r*r*yInit*y*y*y*y*y*y*b-3.0/16.0*b*b*y*r*r*yInit*yInit*yInit*yInit+b*y*r*r*yInit*yInit*yInit*yInit/24.0-7.0/24.0*b*b*b*b*y*r*r*yInit*yInit*yInit*yInit+b*b*b*y*r*r*yInit*yInit*yInit*yInit/3.0-9.0/16.0*b*b*y*r*yInit*yInit+b*y*r*yInit*yInit/8.0+7.0/8.0*b*b*b*y*r*yInit*yInit-b*b*b*b*b*b*y*r*r*yInit*yInit*yInit*yInit/48.0+b*b*b*b*b*y*r*r*yInit*yInit*yInit*yInit/8.0-9.0/16.0*b*b*b*b*y*r*yInit*yInit+b*b*b*b*b*y*r*yInit*yInit/8.0-11.0/16.0*b*b*y-9.0/64.0*b*b*b*b*y+9.0/16.0*b*b*b*y+r*r*r*y*y*y*y*y*b*b*b*b*b*b*b*yInit*yInit/6.0;      
		MapleGenVar14 = 1/yInit/pow(-1.0+b,4.0)/(y*y);      
		double MapleGenVar12 = MapleGenVar13*MapleGenVar14;      
		MapleGenVar13 = -(-396.0*yInit*b*b+144.0*yInit*b-16.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit-96.0*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit-48.0*yInit*yInit*yInit*yInit*yInit*r*r-1248.0*r*r*r*b*b*b*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit+3456.0*r*r*r*b*b*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit-2592.0*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b+192.0*r*r*r*b*b*b*b*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit+768.0*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b+4800.0*r*r*r*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit-5280.0*r*r*r*b*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit+128.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b*b*b*b+128.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b-1120.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b+896.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b-448.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b+896.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b*b-448.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b*b*b-16.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b*b*b*b*b+324.0*yInit*b*b*b-81.0*yInit*b*b*b*b-2520.0*r*r*yInit*yInit*yInit*yInit*yInit*b*b-504.0*b*b*b*b*b*b*r*r*yInit*yInit*yInit*yInit*yInit+2448.0*b*b*b*b*b*r*r*yInit*yInit*yInit*yInit*yInit-1296.0*b*b*b*b*r*yInit*yInit*yInit+288.0*b*b*b*b*b*r*yInit*yInit*yInit+624.0*r*r*yInit*yInit*yInit*yInit*yInit*b-4800.0*r*r*yInit*yInit*yInit*yInit*yInit*b*b*b*b+4800.0*r*r*yInit*yInit*yInit*yInit*yInit*b*b*b-1296.0*r*yInit*yInit*yInit*b*b+288.0*r*yInit*yInit*yInit*b+2016.0*r*b*b*b*yInit*yInit*yInit)/(yInit*yInit*yInit)/pow(-1.0+b,4.0)/1152.0;      
		double MapleGenVar11 = MapleGenVar12+MapleGenVar13;      
		MapleGenVar12 = Time*Time;      
		double MapleGenVar10 = MapleGenVar11*MapleGenVar12;      
		double MapleGenVar8 = MapleGenVar9*MapleGenVar10;      
		double MapleGenVar6 = MapleGenVar7+MapleGenVar8;      
		double MapleGenVar4 = MapleGenVar5+MapleGenVar6;      
		MapleGenVar5 = sqrt(2.0);      
		double MapleGenVar3 = MapleGenVar4*MapleGenVar5;      
		double MapleGenVar1 = MapleGenVar2*MapleGenVar3;      
		MapleGenVar2 = 1/sqrt(3.141592653589793*Time)*exp(-pow(y-yInit,2.0)/Time/2.0)*exp(-1/(-1.0+b)*r*y*y/2.0+1/(-1.0+b)*r*y*y*b-1/(-1.0+b)*r*y*y*b*b/2.0+1/(-1.0+b)*b*log(y)/2.0+1/(-1.0+b)*r*yInit*yInit/2.0-1/(-1.0+b)*r*yInit*yInit*b+1/(-1.0+b)*r*yInit*yInit*b*b/2.0-1/(-1.0+b)*b*log(yInit)/2.0);      
		return MapleGenVar1*MapleGenVar2/Vol/pow(Value,beta);

/*		double MapleGenVar2 = 1.0/2.0;      
		double MapleGenVar5 = 1.0;      
		double MapleGenVar7 = 1/(y-yInit)*(-(24.0*r*r*y*y*y*y*b*b+4.0*r*r*y*y*y*y-16.0*r*r*y*y*y*y*b+4.0*r*r*y*y*y*y*b*b*b*b-16.0*r*r*y*y*y*y*b*b*b+60.0*r*y*y*b*b-48.0*r*y*y*b-24.0*r*b*b*b*y*y+3.0*b*b+12.0*r*y*y-6.0*b)/pow(-1.0+b,2.0)/y/24.0+(24.0*r*r*yInit*yInit*yInit*yInit*b*b+4.0*r*r*yInit*yInit*yInit*yInit-16.0*r*r*yInit*yInit*yInit*yInit*b+4.0*r*r*yInit*yInit*yInit*yInit*b*b*b*b-16.0*r*r*yInit*yInit*yInit*yInit*b*b*b+60.0*r*yInit*yInit*b*b-48.0*r*yInit*yInit*b-24.0*r*b*b*b*yInit*yInit+3.0*b*b+12.0*r*yInit*yInit-6.0*b)/pow(-1.0+b,2.0)/yInit/24.0)*Time;
		double MapleGenVar9 = 1/(pow(y-yInit,2.0));      
		double MapleGenVar15 = -r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit/36.0+r*r*y*y*y*y*y*b*b*b/3.0-7.0/24.0*r*r*y*y*y*y*y*b*b*b*b-yInit*yInit*r*r*y*y*y/12.0-3.0/16.0*r*r*y*y*y*y*y*b*b-r*r*y*y*y*y*y*b*b*b*b*b*b/48.0+r*r*y*y*y*y*y*b*b*b*b*b/8.0-9.0/16.0*b*b*b*b*y*y*y*r-9.0/16.0*r*y*y*y*b*b+r*y*y*y*b/8.0+b*b*b*b*b*y*y*y*r/8.0-r*r*r*y*y*y*y*y*yInit*yInit/12.0+r*r*y*y*y*y*y*b/24.0+7.0/8.0*r*b*b*b*y*y*y-r*r*r*yInit*yInit*yInit*yInit*y*y*y/12.0+11.0/32.0*yInit*b*b-yInit*b/8.0+9.0/128.0*b*b*b*b*yInit-9.0/32.0*b*b*b*yInit+b*y/4.0-9.0/4.0*r*r*r*yInit*yInit*yInit*yInit*b*b*y*y*y;
		double MapleGenVar14 = MapleGenVar15+2.0/3.0*r*r*r*yInit*yInit*yInit*yInit*b*y*y*y+25.0/6.0*r*r*r*b*b*b*yInit*yInit*yInit*yInit*y*y*y-4.0*yInit*yInit*r*r*b*b*y*y*y+yInit*yInit*r*r*b*y*y*y-31.0/4.0*yInit*yInit*r*r*b*b*b*b*y*y*y+23.0/3.0*yInit*yInit*r*r*b*b*b*y*y*y+4.0*r*r*b*b*b*b*b*yInit*yInit*y*y*y-5.0/6.0*r*r*b*b*b*b*b*b*yInit*yInit*y*y*y-55.0/12.0*r*r*r*b*b*b*b*yInit*yInit*yInit*yInit*y*y*y+3.0*r*r*r*b*b*b*b*b*yInit*yInit*yInit*yInit*y*y*y-13.0/12.0*r*r*r*b*b*b*b*b*b*yInit*yInit*yInit*yInit*y*y*y+r*r*r*b*b*b*b*b*b*b*yInit*yInit*yInit*yInit*y*y*y/6.0-35.0/18.0*r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit*b*b*b*b+14.0/9.0*r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit*b*b*b-55.0/12.0*r*r*r*y*y*y*y*y*b*b*b*b*yInit*yInit-9.0/4.0*r*r*r*y*y*y*y*y*yInit*yInit*b*b+2.0/3.0*r*r*r*y*y*y*y*y*yInit*yInit*b+25.0/6.0*r*r*r*y*y*y*y*y*b*b*b*yInit*yInit-7.0/9.0*r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit*b*b+2.0/9.0*r*r*r*r*y*y*y*y*y*yInit*yInit*yInit*yInit*b+14.0/9.0*r*r*r*r*y*y*y*y*y*b*b*b*b*b*yInit*yInit*yInit*yInit+3.0*r*r*r*y*y*y*y*y*b*b*b*b*b*yInit*yInit;
		MapleGenVar15 = MapleGenVar14-7.0/9.0*r*r*r*r*y*y*y*y*y*b*b*b*b*b*b*yInit*yInit*yInit*yInit-13.0/12.0*r*r*r*y*y*y*y*y*b*b*b*b*b*b*yInit*yInit+2.0/9.0*r*r*r*r*y*y*y*y*y*b*b*b*b*b*b*b*yInit*yInit*yInit*yInit+r*r*r*y*y*y*y*y*b*b*b*b*b*b*b*yInit*yInit/6.0-r*r*r*r*y*y*y*y*y*b*b*b*b*b*b*b*b*yInit*yInit*yInit*yInit/36.0+r*r*r*r*yInit*y*y*y*y*y*y*y*y/72.0-7.0/9.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b+35.0/36.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b-r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b*b*b*b/9.0+r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b*b*b*b*b/72.0+7.0/18.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b+7.0/18.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b*b*b-7.0/9.0*r*r*r*r*yInit*y*y*y*y*y*y*y*y*b*b*b*b*b-r*r*r*r*yInit*y*y*y*y*y*y*y*y*b/9.0+yInit*r*r*y*y*y*y/24.0-17.0/8.0*yInit*r*r*y*y*y*y*b*b*b*b*b+7.0/16.0*yInit*r*r*y*y*y*y*b*b*b*b*b*b+35.0/16.0*yInit*r*r*y*y*y*y*b*b-11.0/16.0*b*b*y-9.0/64.0*b*b*b*b*y+9.0/16.0*b*b*b*y;      
		double MapleGenVar13 = MapleGenVar15+25.0/6.0*yInit*r*r*y*y*y*y*b*b*b*b-25.0/6.0*yInit*r*r*y*y*y*y*b*b*b-13.0/24.0*yInit*r*r*y*y*y*y*b+r*r*r*yInit*y*y*y*y*y*y/12.0+55.0/12.0*r*r*r*yInit*y*y*y*y*y*y*b*b*b*b+9.0/4.0*r*r*r*yInit*y*y*y*y*y*y*b*b-25.0/6.0*r*r*r*yInit*y*y*y*y*y*y*b*b*b-3.0*r*r*r*yInit*y*y*y*y*y*y*b*b*b*b*b+13.0/12.0*r*r*r*yInit*y*y*y*y*y*y*b*b*b*b*b*b-r*r*r*yInit*y*y*y*y*y*y*b*b*b*b*b*b*b/6.0-2.0/3.0*r*r*r*yInit*y*y*y*y*y*y*b-3.0/16.0*b*b*y*r*r*yInit*yInit*yInit*yInit+b*y*r*r*yInit*yInit*yInit*yInit/24.0-7.0/24.0*b*b*b*b*y*r*r*yInit*yInit*yInit*yInit+b*b*b*y*r*r*yInit*yInit*yInit*yInit/3.0-9.0/16.0*b*b*y*r*yInit*yInit+b*y*r*yInit*yInit/8.0+7.0/8.0*b*b*b*y*r*yInit*yInit-b*b*b*b*b*b*y*r*r*yInit*yInit*yInit*yInit/48.0+b*b*b*b*b*y*r*r*yInit*yInit*yInit*yInit/8.0-9.0/16.0*b*b*b*b*y*r*yInit*yInit+b*b*b*b*b*y*r*yInit*yInit/8.0;      
		MapleGenVar14 = 1/yInit/pow(-1.0+b,4.0)/(y*y);      
		double MapleGenVar12 = MapleGenVar13*MapleGenVar14;      
		MapleGenVar13 = -(-2520.0*r*r*yInit*yInit*yInit*yInit*yInit*b*b+624.0*r*r*yInit*yInit*yInit*yInit*yInit*b-4800.0*r*r*yInit*yInit*yInit*yInit*yInit*b*b*b*b+4800.0*r*r*yInit*yInit*yInit*yInit*yInit*b*b*b-1296.0*r*yInit*yInit*yInit*b*b+288.0*r*yInit*yInit*yInit*b+2016.0*r*b*b*b*yInit*yInit*yInit-504.0*b*b*b*b*b*b*r*r*yInit*yInit*yInit*yInit*yInit+2448.0*b*b*b*b*b*r*r*yInit*yInit*yInit*yInit*yInit-1296.0*b*b*b*b*r*yInit*yInit*yInit+288.0*b*b*b*b*b*r*yInit*yInit*yInit-396.0*yInit*b*b+144.0*yInit*b-81.0*b*b*b*b*yInit+324.0*b*b*b*yInit-16.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit-48.0*yInit*yInit*yInit*yInit*yInit*r*r-96.0*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit-2592.0*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b+768.0*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b+4800.0*r*r*r*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit-5280.0*r*r*r*b*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit+3456.0*r*r*r*b*b*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit-1248.0*r*r*r*b*b*b*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit+192.0*r*r*r*b*b*b*b*b*b*b*yInit*yInit*yInit*yInit*yInit*yInit*yInit-1120.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b+896.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b-448.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b+128.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b+896.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b*b-448.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b*b*b+128.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b*b*b*b-16.0*r*r*r*r*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*yInit*b*b*b*b*b*b*b*b)/(yInit*yInit*yInit)/pow(-1.0+b,4.0)/1152.0;      
		double MapleGenVar11 = MapleGenVar12+MapleGenVar13;      
		MapleGenVar12 = Time*Time;      
		double MapleGenVar10 = MapleGenVar11*MapleGenVar12;      
		double MapleGenVar8 = MapleGenVar9*MapleGenVar10;      
		double MapleGenVar6 = MapleGenVar7+MapleGenVar8;      
		double MapleGenVar4 = MapleGenVar5+MapleGenVar6;      
		MapleGenVar5 = sqrt(2.0);      
		double MapleGenVar3 = MapleGenVar4*MapleGenVar5;      
		double MapleGenVar1 = MapleGenVar2*MapleGenVar3;      
		MapleGenVar2 = 1/sqrt(0.3141592653589793*Time)*exp(-pow(y-yInit,2.0)/Time/2.0)*exp(-1/(-1.0+b)*r*y*y/2.0+1/(-1.0+b)*r*y*y*b-1/(-1.0+b)*r*y*y*b*b/2.0+1/(-1.0+b)*b*log(y)/2.0+1/(-1.0+b)*r*yInit*yInit/2.0-1/(-1.0+b)*r*yInit*yInit*b+1/(-1.0+b)*r*yInit*yInit*b*b/2.0-1/(-1.0+b)*b*log(yInit)/2.0);      
		double Result=  MapleGenVar1*MapleGenVar2/Vol/pow(Value,beta);
		return Result; */
	}
}


#endif
