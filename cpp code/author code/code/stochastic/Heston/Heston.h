#ifndef HESTON_H
#define HESTON_H

#include "SVProcess.h" 
#include "MiscMaths.h"

class HestonProcess : public SVProcess
{
public:
	HestonProcess(double LogAsset_, double Var_, double rate_, double kappa_, double theta_, double sigma_, double corr_);
	virtual void Evolve(double Time, MTRand& MT); 
	virtual void Evolve(double Time, unsigned long N, MTRand& MT);
	virtual double InitsLow(double Time, double Deviations) const;
	virtual double InitsHigh(double Time, double Deviations) const;
	virtual double InitvLow(double Time, double Deviations) const;
	virtual double InitvHigh(double Time, double Deviations) const;
	virtual double InitsLow(double sValue, double Time, double Deviations) const;
	virtual double InitsHigh(double sValue, double Time, double Deviations) const;
	virtual double InitvLow(double vValue, double Time, double Deviations) const;
	virtual double InitvHigh(double vValue, double Time, double Deviations) const;
	virtual double Density(double Time, double AssetValue, double VarValue) const;
	virtual double Density(double Time, double AssetValue, double VarValue, double AssetBase, double VarBase) const;
	virtual double Density1(double Time, double AssetValue, double VarValue, double AssetBase, double VarBase) const;
	virtual double Density2(double Time, double AssetValue, double VarValue, double AssetBase, double VarBase) const;
	virtual double Density3(double Time, double AssetValue, double VarValue, double AssetBase, double VarBase) const;
	virtual double GetValue() const;
	virtual double GetRate() const;
	virtual double GetVar() const;
	virtual void SetValues(double NewLogAsset, double NewVar); 
	bool FellerOK() const;
private:
	double LogAsset;
	double Var;
	double rate;
	double kappa;
	double theta;
	double sigma;
	double corr;
}; 

HestonProcess::HestonProcess(double LogAsset_, double Var_, double rate_, double kappa_, double theta_, double sigma_, double corr_) : LogAsset(LogAsset_),
Var(Var_), kappa(kappa_), theta(theta_), sigma(sigma_), rate(rate_), corr(corr_)
{
}

void HestonProcess::Evolve(double Time, MTRand& MT)
{
	double r1 = MT.randNorm();
	LogAsset += (rate - 0.5*Var)*Time + sqrt(Var*Time)*r1;
	if (Var > 0)
	{
		double r2 = sqrt(1-corr*corr)*MT.randNorm() + corr*r1; 
		Var += kappa*(theta-Var)*Time + sigma*sqrt(Var*Time)*r2; 
	}
	if (Var < 0)
	{
		Var = 0;
	}
}

void HestonProcess::Evolve(double Time, unsigned long N, MTRand& MT)
{
	double TimeStep = Time/static_cast<double>(N);
	for (unsigned long i = 0; i < N; i++)
	{
		Evolve(TimeStep,MT);
	}
}

double HestonProcess::GetRate() const
{
	return rate;
}

double HestonProcess::GetValue() const
{
	return LogAsset; 
}

double HestonProcess::GetVar() const
{
	return Var;
}

void HestonProcess::SetValues(double NewLogAsset, double NewVar)
{
	LogAsset = NewLogAsset;
	Var = NewVar;
}

double HestonProcess::Density(double Time, double AssetValue, double VarValue) const
{
	return Density(Time,AssetValue,VarValue,LogAsset,Var);
}

double HestonProcess::Density(double Time, double AssetValue, double VarValue, double AssetBase, double VarBase) const
{
	return Density3(Time,AssetValue,VarValue,AssetBase,VarBase);
}

double HestonProcess::Density1(double Time, double AssetValue, double VarValue, double AssetBase, double VarBase) const
{
	double m = rate;
	double s = -sigma;
	double r = -sqrt(1-corr*corr);
	double a = kappa*theta;
	double b = kappa;
	double X1 = AssetValue-AssetBase;
	double X2 = VarValue-VarBase;

	double MapleGenVar3 = -log(2.0*Time*3.14159265358979323846)-log(r*s*VarValue);
	double MapleGenVar4 = MapleGenVar3;
	double MapleGenVar8 = -7.0/11520.0*s*s*s*s*X1*X1*X1*X1*X1*X1/(r*r*r*r*r*r)/(VarBase*VarBase*VarBase*VarBase*VarBase)+7.0/1920.0*sqrt(1.0-r*r)*s*s*s*X1*X1*X1*X1*X1*X2/(r*r*r*r*r*r)/(VarBase*VarBase*VarBase*VarBase*VarBase)+(-105.0+298.0*r*r)*s*s*X1*X1*X1*X1*X2*X2/(r*r*r*r*r*r)/(VarBase*VarBase*VarBase*VarBase*VarBase)/11520.0+(35.0-228.0*r*r)*sqrt(1.0-r*r)*s*X1*X1*X1*X2*X2*X2/(r*r*r*r*r*r)/(VarBase*VarBase*VarBase*VarBase*VarBase)/2880.0-(105.0-1368.0*r*r+2008.0*r*r*r*r)*X1*X1*X2*X2*X2*X2/(r*r*r*r*r*r)/(VarBase*VarBase*VarBase*VarBase*VarBase)/11520.0+sqrt(1.0-r*r)*(21.0-428.0*r*r+1152.0*r*r*r*r)*X1*X2*X2*X2*X2*X2/(r*r*r*r*r*r)/s/(VarBase*VarBase*VarBase*VarBase*VarBase)/5760.0-(7.0-214.0*r*r+1152.0*r*r*r*r)*X2*X2*X2*X2*X2*X2/(r*r*r*r*r*r)/(s*s)/(VarBase*VarBase*VarBase*VarBase*VarBase)/11520.0-s*s*X1*X1*X1*X1*X2/(r*r*r*r)/(VarBase*VarBase*VarBase*VarBase)/64.0+sqrt(1.0-r*r)*s*X1*X1*X1*X2*X2/(r*r*r*r)/(VarBase*VarBase*VarBase*VarBase)/16.0+(-3.0+6.0*r*r)*X1*X1*X2*X2*X2/(r*r*r*r)/(VarBase*VarBase*VarBase*VarBase)/32.0;
	double MapleGenVar7 = MapleGenVar8+(1.0-4.0*r*r)*sqrt(1.0-r*r)*X1*X2*X2*X2*X2/(r*r*r*r)/s/(VarBase*VarBase*VarBase*VarBase)/16.0+(-1.0+8.0*r*r)*X2*X2*X2*X2*X2/(r*r*r*r)/(s*s)/(VarBase*VarBase*VarBase*VarBase)/64.0+s*s*X1*X1*X1*X1/(r*r*r*r)/(VarBase*VarBase*VarBase)/96.0-sqrt(1.0-r*r)*s*X1*X1*X1*X2/(r*r*r*r)/(VarBase*VarBase*VarBase)/24.0+(3.0-10.0*r*r)*X1*X1*X2*X2/(r*r*r*r)/(VarBase*VarBase*VarBase)/48.0+sqrt(1.0-r*r)*(-1.0+8.0*r*r)*X1*X2*X2*X2/(r*r*r*r)/s/(VarBase*VarBase*VarBase)/24.0+(1.0-16.0*r*r)*X2*X2*X2*X2/(r*r*r*r)/(s*s)/(VarBase*VarBase*VarBase)/96.0+X1*X1*X2/(r*r)/(VarBase*VarBase)/4.0-sqrt(1.0-r*r)*X1*X2*X2/(r*r)/s/(VarBase*VarBase)/2.0+X2*X2*X2/(r*r)/(s*s)/(VarBase*VarBase)/4.0-(s*s*X1*X1-2.0*sqrt(1.0-r*r)*s*X1*X2+X2*X2)/(r*r)/(s*s)/VarBase/2.0;      
	MapleGenVar8 = 1/Time;      
	double MapleGenVar6 = MapleGenVar7*MapleGenVar8;      
	MapleGenVar7 = 7.0/1920.0*s*s*s*s*X1*X1*X1*X1/(r*r*r*r)/(VarBase*VarBase*VarBase*VarBase);      
	double MapleGenVar5 = MapleGenVar6+MapleGenVar7;      
	double MapleGenVar2 = MapleGenVar4+MapleGenVar5;      
	double MapleGenVar1 = MapleGenVar2-1/(r*r*r*r)/(VarBase*VarBase*VarBase*VarBase)*s*(30.0*a*sqrt(1.0-r*r)+s*(-30.0*m+7.0*sqrt(1.0-r*r)*s))*X1*X1*X1*X2/480.0+1/(r*r*r*r)/(VarBase*VarBase*VarBase*VarBase)*(-540.0*a*(-1.0+r*r)+s*(-540.0*m*sqrt(1.0-r*r)+(63.0-160.0*r*r)*s))*X1*X1*X2*X2/2880.0+1/(r*r*r*r)/s/(VarBase*VarBase*VarBase*VarBase)*(270.0*a*sqrt(1.0-r*r)*(-1.0+2.0*r*r)+s*(-270.0*m*(-1.0+2.0*r*r)+sqrt(1.0-r*r)*(-21.0+118.0*r*r)*s))*X1*X2*X2*X2/1440.0+1/(r*r*r*r)/(s*s)/(VarBase*VarBase*VarBase*VarBase)*(-360.0*a*(-1.0+5.0*r*r)+s*(360.0*m*sqrt(1.0-r*r)*(-1.0+4.0*r*r)+(21.0-236.0*r*r)*s))*X2*X2*X2*X2/5760.0+s*(a*sqrt(1.0-r*r)-m*s)*X1*X1*X1/(r*r*r*r)/(VarBase*VarBase*VarBase)/24.0;      
	MapleGenVar2 = MapleGenVar1+(3.0*a*(-1.0+r*r)+s*(3.0*m*sqrt(1.0-r*r)+r*r*s))*X1*X1*X2/(r*r*r*r)/(VarBase*VarBase*VarBase)/24.0+1/(r*r*r*r)/s/(VarBase*VarBase*VarBase)*(a*(3.0-10.0*r*r)*sqrt(1.0-r*r)+s*(m*(-3.0+10.0*r*r)-2.0*r*r*sqrt(1.0-r*r)*s))*X1*X2*X2/24.0+(a*(-1.0+9.0*r*r)+s*(m*(1.0-8.0*r*r)*sqrt(1.0-r*r)+r*r*s))*X2*X2*X2/(r*r*r*r)/(s*s)/(VarBase*VarBase*VarBase)/24.0-s*s*X1*X1/(r*r)/(VarBase*VarBase)/24.0;      
	MapleGenVar3 = MapleGenVar2+(6.0*a*sqrt(1.0-r*r)-6.0*m*s+sqrt(1.0-r*r)*s*s)*X1*X2/(r*r)/s/(VarBase*VarBase)/12.0-(12.0*a-12.0*m*sqrt(1.0-r*r)*s+s*s)*X2*X2/(r*r)/(s*s)/(VarBase*VarBase)/24.0;      
	MapleGenVar4 = MapleGenVar3-X1*(2.0*a*sqrt(1.0-r*r)-2.0*m*s-2.0*b*sqrt(1.0-r*r)*VarBase+s*VarBase)/(r*r)/s/VarBase/2.0;      
	MapleGenVar5 = MapleGenVar4+X2*(2.0*a-2.0*m*sqrt(1.0-r*r)*s-2.0*b*VarBase+sqrt(1.0-r*r)*s*VarBase)/(r*r)/(s*s)/VarBase/2.0;      
	MapleGenVar6 = MapleGenVar5;      
	double MapleGenVar9 = Time;      
	double MapleGenVar11 = s*(-a*sqrt(1.0-r*r)+m*s)*X1/(r*r)/(VarBase*VarBase)/12.0+1/(r*r*r*r)/(VarBase*VarBase*VarBase)*X1*X1*(-60.0*a*a*(-3.0+2.0*r*r)+180.0*m*m*s*s+2.0*r*r*s*s*s*s-60.0*a*s*(6.0*m*sqrt(1.0-r*r)+r*r*s)-60.0*b*b*VarBase*VarBase+60.0*b*sqrt(1.0-r*r)*s*VarBase*VarBase-15.0*s*s*VarBase*VarBase)/2880.0+1/(r*r)/(s*s)/(VarBase*VarBase)*X2*(12.0*a*a+12.0*m*m*s*s-4.0*m*sqrt(1.0-r*r)*s*s*s+2.0*r*r*s*s*s*s-4.0*a*s*(6.0*m*sqrt(1.0-r*r)+(-1.0+4.0*r*r)*s)-12.0*b*b*VarBase*VarBase+12.0*b*sqrt(1.0-r*r)*s*VarBase*VarBase-3.0*s*s*VarBase*VarBase)/48.0;      
	double MapleGenVar12 = MapleGenVar11+1/(r*r*r*r)/(s*s)/(VarBase*VarBase*VarBase)*X2*X2*(-180.0*a*a*(-1.0+4.0*r*r)-60.0*m*m*(-3.0+10.0*r*r)*s*s+240.0*m*r*r*sqrt(1.0-r*r)*s*s*s+2.0*r*r*s*s*s*s-96.0*r*r*r*r*s*s*s*s+60.0*a*s*(2.0*m*sqrt(1.0-r*r)*(-3.0+10.0*r*r)+r*r*(-5.0+14.0*r*r)*s)-60.0*b*b*VarBase*VarBase+120.0*b*b*r*r*VarBase*VarBase+60.0*b*sqrt(1.0-r*r)*s*VarBase*VarBase-120.0*b*r*r*sqrt(1.0-r*r)*s*VarBase*VarBase-15.0*s*s*VarBase*VarBase+30.0*r*r*s*s*VarBase*VarBase)/2880.0;      
	double MapleGenVar10 = MapleGenVar12+1/(r*r*r*r)/s/(VarBase*VarBase*VarBase)*X1*X2*(60.0*a*a*sqrt(1.0-r*r)*(-3.0+2.0*r*r)-180.0*m*m*sqrt(1.0-r*r)*s*s-120.0*m*r*r*s*s*s-2.0*r*r*sqrt(1.0-r*r)*s*s*s*s+180.0*a*s*(-2.0*m*(-1.0+r*r)+r*r*sqrt(1.0-r*r)*s)+60.0*b*b*sqrt(1.0-r*r)*VarBase*VarBase-60.0*b*s*VarBase*VarBase+60.0*b*r*r*s*VarBase*VarBase+15.0*sqrt(1.0-r*r)*s*s*VarBase*VarBase)/1440.0-1/(r*r)/(s*s)/VarBase*(12.0*a*a+12.0*m*m*s*s+2.0*r*r*s*s*s*s-12.0*m*s*(-2.0*b*sqrt(1.0-r*r)+s)*VarBase+12.0*b*b*VarBase*VarBase-12.0*b*sqrt(1.0-r*r)*s*VarBase*VarBase+3.0*s*s*VarBase*VarBase-12.0*a*(2.0*m*sqrt(1.0-r*r)*s+r*r*s*s+2.0*b*VarBase-sqrt(1.0-r*r)*s*VarBase))/24.0;      
	MapleGenVar8 = MapleGenVar9*MapleGenVar10;      
	MapleGenVar9 = -Time*Time*(60.0*a*a*(1.0+r*r)+60.0*m*m*s*s+23.0*r*r*s*s*s*s-120.0*a*s*(m*sqrt(1.0-r*r)+r*r*s))/(r*r)/(VarBase*VarBase)/1440.0;      
	MapleGenVar7 = MapleGenVar8+MapleGenVar9;
	return exp(MapleGenVar6+MapleGenVar7);
}

double HestonProcess::Density2(double Time, double AssetValue, double VarValue, double AssetBase, double VarBase) const
{
	double m = rate;
	double s = -sigma;
	double r = -sqrt(1-corr*corr);
	double a = kappa*theta;
	double b = kappa;
	double X1 = AssetValue-AssetBase;
	double X2 = VarValue-VarBase;

	double r2 = r*r;
	double vI2 = VarBase*VarBase;
	double vI3 = vI2*VarBase;
	double s2 = s*s;
	double X12 = X1*X1;
	double X22 = X2*X2;

	double CMinus = -((7*s2*s2*X12*X12*X12)/(11520*r2*r2*r2*vI2*vI3)) + (7*sqrt(1 - r2)*s2*s*X12*X12*X1*X2)/(1920*r2*r2*r2*vI2*vI3) + ((-105 + 298*r2)*s2*X12*X12*X22)/(11520*r2*r2*r2*vI2*vI3) + ((35 - 228*r2)*sqrt(1 - r2)*s*X12*X1*X22*X2)/(2880*r2*r2*r2*vI2*vI3) - ((105 - 1368*r2 + 2008*r2*r2)*X12*X22*X22)/(11520*r2*r2*r2*vI2*vI3) + (sqrt(1 - r2)*(21 - 428*r2 + 1152*r2*r2)*X1*X22*X22*X2)/(5760*r2*r2*r2*s*vI2*vI3) - ((7 - 214*r2 + 1152*r2*r2)*X22*X22*X22)/(11520*r2*r2*r2*s2*vI2*vI3) - (s2*X12*X12*X2)/(64*r2*r2*vI2*vI2) + (sqrt(1 - r2)*s*X12*X1*X22)/(16*r2*r2*vI2*vI2) + (3*(-1 + 2*r2)*X12*X22*X2)/(32*r2*r2*vI2*vI2) + ((1 - 4*r2)*sqrt(1 - r2)*X1*X22*X22)/(16*r2*r2*s*vI2*vI2) + ((-1 + 8*r2)*X22*X22*X2)/(64*r2*r2*s2*vI2*vI2) + (s2*X12*X12)/(96*r2*r2*vI3) - (sqrt(1 - r2)*s*X12*X1*X2)/(24*r2*r2*vI3) + ((3 - 10*r2)*X12*X22)/(48*r2*r2*vI3) + (sqrt(1 - r2)*(-1 + 8*r2)*X1*X22*X2)/(24*r2*r2*s*vI3) + ((1 - 16*r2)*X22*X22)/(96*r2*r2*s2*vI3) + (X12*X2)/(4*r2*vI2) - (sqrt(1 - r2)*X1*X22)/(2*r2*s*vI2) + X22*X2/(4*r2*s2*vI2) - (s2*X12 - 2*sqrt(1 - r2)*s*X1*X2 + X22)/(2*r2*s2*VarBase);
	double CZero = (7*s2*s2*X12*X12)/(1920*r2*r2*vI2*vI2)-(1/(480*r2*r2*vI2*vI2))*(s*(30*a*sqrt(1-r2)+s*(-30*m+7*sqrt(1-r2)*s))*X12*X1*X2)+(1/(2880*r2*r2*vI2*vI2))*((-540*a*(-1+r2)+s*(-540*m*sqrt(1-r2)+(63-160*r2)*s))*X12*X22)+(1/(1440*r2*r2*s*vI2*vI2))*((270*a*sqrt(1-r2)*(-1+2*r2)+s*(-270*m*(-1+2*r2)+sqrt(1-r2)*(-21+118*r2)*s))*X1*X22*X2)+(1/(5760*r2*r2*s2*vI2*vI2))*((-360*a*(-1+5*r2)+s*(360*m*sqrt(1-r2)*(-1+4*r2)+(21-236*r2)*s))*X22*X22)+(s*(a*sqrt(1-r2)-m*s)*X12*X1)/(24*r2*r2*vI3)+((3*a*(-1+r2)+s*(3*m*sqrt(1-r2)+r2*s))*X12*X2)/(24*r2*r2*vI3)+(1/(24*r2*r2*s*vI3))*((a*(3-10*r2)*sqrt(1-r2)+s*(m*(-3+10*r2)-2*r2*sqrt(1-r2)*s))*X1*X22)+((a*(-1+9*r2)+s*(m*(1-8*r2)*sqrt(1-r2)+r2*s))*X22*X2)/(24*r2*r2*s2*vI3)-(s2*X12)/(24*r2*vI2)+((6*a*sqrt(1-r2)-6*m*s+sqrt(1-r2)*s2)*X1*X2)/(12*r2*s*vI2)-((12*a-12*m*sqrt(1-r2)*s+s2)*X22)/(24*r2*s2*vI2)-(X1*(2*a*sqrt(1-r2)-2*m*s-2*b*sqrt(1-r2)*VarBase+s*VarBase))/(2*r2*s*VarBase)+(X2*(2*a-2*m*sqrt(1-r2)*s-2*b*VarBase+sqrt(1-r2)*s*VarBase))/(2*r2*s2*VarBase);
	double COne = (s*((-a)*sqrt(1-r2)+m*s)*X1)/(12*r2*vI2)+(1/(2880*r2*r2*vI3))*(X12*(-60*a*a*(-3+2*r2)+180*m*m*s2+2*r2*s2*s2-60*a*s*(6*m*sqrt(1-r2)+r2*s)-60*b*b*vI2+60*b*sqrt(1-r2)*s*vI2-15*s2*vI2))+(1/(48*r2*s2*vI2))*(X2*(12*a*a+12*m*m*s2-4*m*sqrt(1-r2)*s2*s+2*r2*s2*s2-4*a*s*(6*m*sqrt(1-r2)+(-1+4*r2)*s)-12*b*b*vI2+12*b*sqrt(1-r2)*s*vI2-3*s2*vI2))+(1/(2880*r2*r2*s2*vI3))*(X22*(-180*a*a*(-1+4*r2)-60*m*m*(-3+10*r2)*s2+240*m*r2*sqrt(1-r2)*s2*s+2*r2*s2*s2-96*r2*r2*s2*s2+60*a*s*(2*m*sqrt(1-r2)*(-3+10*r2)+r2*(-5+14*r2)*s)-60*b*b*vI2+120*b*b*r2*vI2+60*b*sqrt(1-r2)*s*vI2-120*b*r2*sqrt(1-r2)*s*vI2-15*s2*vI2+30*r2*s2*vI2))+(1/(1440*r2*r2*s*vI3))*(X1*X2*(60*a*a*sqrt(1-r2)*(-3+2*r2)-180*m*m*sqrt(1-r2)*s2-120*m*r2*s2*s-2*r2*sqrt(1-r2)*s2*s2+180*a*s*(-2*m*(-1+r2)+r2*sqrt(1-r2)*s)+60*b*b*sqrt(1-r2)*vI2-60*b*s*vI2+60*b*r2*s*vI2+15*sqrt(1-r2)*s2*vI2))-(1/(24*r2*s2*VarBase))*(12*a*a+12*m*m*s2+2*r2*s2*s2-12*m*s*(-2*b*sqrt(1-r2)+s)*VarBase+12*b*b*vI2-12*b*sqrt(1-r2)*s*vI2+3*s2*vI2-12*a*(2*m*sqrt(1-r2)*s+r2*s2+2*b*VarBase-sqrt(1-r2)*s*VarBase));
	double CTwo = -((60*a*a*(1+r2)+60*m*m*s2+23*r2*s2*s2-120*a*s*(m*sqrt(1-r2)+r2*s))/(720*r2*vI2));
	return exp(CMinus/Time + CZero + Time*COne + 0.5*Time*Time*CTwo)/r/s/VarValue/2.0/3.14159265358979323846/Time;
}


double HestonProcess::Density3(double Time, double AssetValue, double VarValue, double AssetBase, double VarBase) const
{
	double x1 = AssetValue;
	double x2 = VarValue;
	double x10 = AssetBase;
	double x20 = VarBase;
	double a1 = rate;
	double a2 = kappa*theta;
	double b2 = -kappa;
	double morr = 1-corr*corr;

	double CMinus = ((pow(x2,2)*(-15+16*pow(corr,2))+pow(x20,2)*(-15-72*morr+16*pow(corr,2))+2*x2*(-x1+x10)*corr*sigma+pow(x1-x10,2)*pow(sigma,2)+2*x20*(x2*(15+12*morr-16*pow(corr,2))+(x1-x10)*corr*sigma))*(pow(x20,2)+pow(x2,2)+2*x2*(-x1+x10)*corr*sigma+pow(x1-x10,2)*pow(sigma,2)-2*x20*(x2 + (-x1 + x10)*corr*sigma)))/(96*pow(morr,2)*pow(x20,3)*pow(sigma,2));
	double CZero = -(pow(x1 - x10,2)*pow(sigma,4)+24*x20*(x1-x10)*sigma*(a2*corr+b2*x20*corr-a1*sigma+0.5*x20*sigma)+24*x20*(x20-x2)*(a2+b2*x20 +(-a1+0.5*x20)*corr*sigma)+pow(x20-x2,2)*(12*a2+sigma*(-12*a1*corr+sigma))+2*(x20-x2)*(x1-x10)*sigma*(6*a2*corr+sigma*(-6*a1+corr*sigma)))/(24*morr*pow(x20,2)*pow(sigma,2));
	double COne = (-0.5*pow(a2,2) - 0.5*pow(b2,2)*pow(x20,2) + b2*(a1 - 0.5*x20)*x20*corr*sigma+pow(sigma,2)*(-0.5*pow(a1,2)+0.5*a1*x20-0.125*pow(x20,2)+((-1 + pow(corr,2))/12)*pow(sigma,2))+a2*(-1*b2*x20+sigma*(a1*corr-0.5*x20*corr+0.5*sigma-0.5*pow(corr,2)*sigma)))/(morr*x20*pow(sigma,2));
	double CTwo = -(60*pow(a2,2)*(1 + pow(corr,2)) + 60*pow(a1,2)*pow(sigma,2)+23*pow(corr,2)*pow(sigma,4)-120*a2*sigma*(a1*sqrt(morr)+pow(corr,2)*sigma))/(720*pow(x20,2)*pow(corr,2));
	return exp(CMinus/Time + CZero + Time*COne + 0.5*Time*Time*CTwo)/sqrt(morr)/sigma/VarValue/2.0/3.14159265358979323846/Time;
}
bool HestonProcess::FellerOK() const
{
	double Proxy = 2.0*kappa*theta/sigma/sigma;
	if (Proxy < 1.0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

double HestonProcess::InitsLow(double Time, double Deviations) const
{
	return LogAsset + (rate - 0.5*Var)*Time - Deviations*sqrt(Var)*sqrt(Time);
}

double HestonProcess::InitsHigh(double Time, double Deviations) const
{
	return LogAsset + (rate - 0.5*Var)*Time + Deviations*sqrt(Var)*sqrt(Time);
}


double HestonProcess::InitsLow(double sValue, double Time, double Deviations) const
{
	return sValue + (rate - 0.5*Var)*Time - Deviations*sqrt(Var)*sqrt(Time);
}

double HestonProcess::InitsHigh(double sValue, double Time, double Deviations) const
{
	return sValue + (rate - 0.5*Var)*Time + Deviations*sqrt(Var)*sqrt(Time);
}


double HestonProcess::InitvLow(double Time, double Deviations) const
{
	return InitvLow(Var,Time,Deviations);
}

double HestonProcess::InitvHigh(double Time, double Deviations) const
{
	return InitvHigh(Var,Time,Deviations);
}

double HestonProcess::InitvLow(double vValue, double Time, double Deviations) const
{
	double Guess = vValue - Deviations*sqrt(vValue*Time);
	return max(Guess,0.091);
}

double HestonProcess::InitvHigh(double vValue, double Time, double Deviations) const
{
	double Guess = vValue + Deviations*sqrt(vValue*Time);
	return min(Guess,0.48);
}


#endif
