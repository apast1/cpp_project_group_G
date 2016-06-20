#ifndef SABR_H
#define SABR_H

#include "SVProcess.h" 
#include "MiscMaths.h"

class SABRProcess : public SVProcess
{
public:
	SABRProcess(double Asset_, double Vol_, double r_, double b_, double a_, double p_);
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
//	virtual double Density1(double Time, double AssetValue, double VolValue, double AssetBase, double VolBase) const;
	virtual double Density2(double Time, double AssetValue, double VolValue, double AssetBase, double VolBase) const;
	virtual double GetValue() const;
	virtual double GetRate() const;
	virtual double GetVar() const;
	virtual void SetValues(double NewLogAsset, double NewVol); 
private:
	double LogAsset;
	double Vol;
	double r;
	double b;
	double a;
	double p;
}; 

SABRProcess::SABRProcess(double LogAsset_, double Vol_, double r_, double b_, double a_, double p_) : LogAsset(LogAsset_), Vol(Vol_), 
r(r_), b(b_), a(a_), p(p_)
{
}


double SABRProcess::GetRate() const
{
	return r;
}

double SABRProcess::GetValue() const
{
	return LogAsset; 
}

double SABRProcess::GetVar() const
{
	return Vol;
}

void SABRProcess::SetValues(double NewLogAsset, double NewVol)
{
	LogAsset = NewLogAsset;
	Vol = NewVol;
}

double SABRProcess::Density(double Time, double AssetValue, double VolValue) const
{
	return Density(Time,AssetValue,VolValue,LogAsset,Vol);
}

double SABRProcess::Density(double Time, double AssetValue, double VolValue, double AssetBase, double VolBase) const
{
	return Density2(Time,AssetValue,VolValue,AssetBase,VolBase);
}

/*
double SABRProcess::Density1(double Time, double AssetValue, double VolValue, double AssetBase, double VolBase) const
{
//  HL density in c form
}
*/

double SABRProcess::Density2(double Time, double AssetValue, double VolValue, double AssetBase, double VolBase) const
{
//  AS density in c form 

	double x1 = AssetValue;
	double x10 = AssetBase;
	double x2 = VolValue;
	double x20 = VolBase;

///////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////


	double cMinusOne = -((p*(16-45*p*p+30*p*p*p*p)*(x1-x10)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20))/(pow(x10,b)*(30*a*(-1+p*p)*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20*x20)))+((137-315*p*p+180*p*p*p*p)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20))/(360*a*a*(-1+p*p)*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20*x20)+(p*(-2+3*p*p)*(x1-x10)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20))/(pow(x10,b)*(3*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20))+((5-6*p*p)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20))/(12*a*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20)-(p*(-5+6*p*p)*(x1-x10)*(x2-x20)*(x2-x20)*(x2-x20))/(pow(x10,b)*(6*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20))+((-11+12*p*p)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20))/(24*a*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20)+(p*(x1-x10)*(x2-x20)*(x2-x20))/(pow(x10,b)*(a*(-1+p*p)*x20*x20*x20))-(x2-x20)*(x2-x20)*(x2-x20)/(2*a*a*(-1+p*p)*x20*x20*x20)-(b*(x1-x10)*(x1-x10)*(x1-x10)*pow(x10,(-1-2*b)))/(2*(-1+p*p)*x20*x20)+(a*a*(x1-x10)*(x1-x10)-2*a*p*(x1-x10)*pow(x10,b)*(x2-x20)+pow(x10,2*b)*(x2-x20)*(x2-x20))/(pow(x10,2*b)*(2*a*a*(-1+p*p)*x20*x20))+((x1-x10)*(x1-x10)*pow(x10,-1-2*b)*(x2-x20)*((-a)*x10+b*p*pow(x10,b)*x20))/(2*a*(-1+p*p)*x20*x20*x20)+((x1-x10)*(x1-x10)*pow(x10,-1-2*b)*(x2-x20)*(x2-x20)*(a*(-5+8*p*p)*x10-6*b*p*(-1+p*p)*pow(x10,b)*x20))/(12*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20)+((x1-x10)*(x1-x10)*pow(x10,-1-2*b)*(x2-x20)*(x2-x20)*(x2-x20)*(2*a*(2-5*p*p)*x10+b*p*(-5+6*p*p)*pow(x10,b)*x20))/(12*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20)+((x1-x10)*(x1-x10)*pow(x10,-1-2*b)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20)*(a*(16-71*p*p+60*p*p*p*p)*x10-10*b*p*(2-5*p*p+3*p*p*p*p)*pow(x10,b)*x20))/(60*a*(-1+p*p)*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20*x20)+((x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-2-4*b)*(a*a*x10*x10+b*(4+7*b)*(-1+p*p)*pow(x10,2*b)*x20*x20))/(24*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20)-(b*(x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-3-4*b)*(2*a*a*x10*x10+(2+5*b+3*b*b)*(-1+p*p)*pow(x10,2*b)*x20*x20))/(24*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20)+((x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-2-3*b)*(x2-x20)*((-a*a)*p*x10*x10+3*a*b*(-1+p*p)*pow(x10,1+b)*x20-b*(1+b)*p*(-1+p*p)*pow(x10,2*b)*x20*x20))/(6*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20)+((x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-2-3*b)*(x2-x20)*(x2-x20)*(4*a*a*p*x10*x10-a*b*(-5+8*p*p)*pow(x10,1+b)*x20+2*b*(1+b)*p*(-1+p*p)*pow(x10,2*b)*x20*x20))/(12*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20)+((x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-2-3*b)*(x2-x20)*(x2-x20)*(x2-x20)*(2*a*a*p*(39-49*p*p)*x10*x10+30*a*b*(2-7*p*p+5*p*p*p*p)*pow(x10,1+b)*x20-5*b*(1+b)*p*(5-11*p*p+6*p*p*p*p)*pow(x10,2*b)*x20*x20))/(180*a*(-1+p*p)*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20*x20)+((x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-3-4*b)*(x2-x20)*(-2*a*a*a*x10*x10*x10+6*a*a*b*p*pow(x10,2+b)*x20-a*b*(4+7*b)*(-1+p*p)*pow(x10,1+2*b)*x20*x20+b*(2+3*b+b*b)*p*(-1+p*p)*pow(x10,3*b)*x20*x20*x20))/(24*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20)+(1/(720*a*(-1+p*p)*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20*x20))*((x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-3-4*b)*(x2-x20)*(x2-x20)*(6*a*a*a*(-13+23*p*p)*x10*x10*x10-360*a*a*b*p*(-1+p*p)*pow(x10,2+b)*x20+5*a*b*(4+7*b)*(5-13*p*p+8*p*p*p*p)*pow(x10,1+2*b)*x20*x20-30*b*(2+3*b+b*b)*p*(-1+p*p)*(-1+p*p)*pow(x10,3*b)*x20*x20*x20))+((x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-4-6*b)*(4*a*a*a*a*x10*x10*x10*x10+5*a*a*b*(4+13*b)*(-1+p*p)*pow(x10,2+2*b)*x20*x20+b*(36+106*b+101*b*b+31*b*b*b)*(-1+p*p)*(-1+p*p)*pow(x10,4*b)*x20*x20*x20*x20))/(720*(-1+p*p)*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20*x20)-(1/(120*a*(-1+p*p)*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20*x20))*((x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-4-5*b)*(x2-x20)*(4*a*a*a*a*p*x10*x10*x10*x10-20*a*a*a*b*(-1+p*p)*pow(x10,3+b)*x20+5*a*a*b*(2+5*b)*p*(-1+p*p)*pow(x10,2+2*b)*x20*x20-5*a*b*(2+5*b+3*b*b)*(-1+p*p)*(-1+p*p)*pow(x10,1+3*b)*x20*x20*x20+b*(6+11*b+6*b*b+b*b*b)*p*(-1+p*p)*(-1+p*p)*pow(x10,4*b)*x20*x20*x20*x20));

	double cZero=-(((x1-x10)*(2*r*pow(x10,1-2*b)-(b*x20*x20)/x10))/(2*(-1+p*p)*x20*x20))-((x2-x20)*(2*p*r*pow(x10,1-b)-b*p*pow(x10,-1+b)*x20*x20))/(2*a*x20*x20-2*a*p*p*x20*x20)+((x1-x10)*(x1-x10)*(a*a*x10*x10+6*(-1+2*b)*r*x10*x10-3*b*pow(x10,2*b)*x20*x20))/(pow(x10,2*(1+b))*(12*(-1+p*p)*x20*x20))+((x2-x20)*(x2-x20)*(12*p*r*pow(x10,1-b)-a*x20))/(12*a*x20*x20*x20-12*a*p*p*x20*x20*x20)+((x1-x10)*(x2-x20)*(12*a*r*x10*x10*x10-2*a*a*p*pow(x10,2+b)*x20-3*(-1+b)*p*pow(x10,b)*x20*(2*r*x10*x10+b*pow(x10,2*b)*x20*x20)))/(pow(x10,2*(1+b))*(12*a*(-1+p*p)*x20*x20*x20))+(pow(x10,-2-b)*(x2-x20)*(x2-x20)*(x2-x20)*(4*a*p*(-5+6*p*p)*r*x10*x10*x10-2*a*a*(-1+p*p)*pow(x10,2+b)*x20+(-1+b)*p*p*pow(x10,b)*x20*(2*r*x10*x10+b*pow(x10,2*b)*x20*x20)))/(24*a*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20)+(pow(x10,-2-b)*(x2-x20)*(x2-x20)*(x2-x20)*(x2-x20)*(-240*a*p*(-2+3*p*p)*r*x10*x10*x10+a*a*(-53+60*p*p)*pow(x10,2+b)*x20-15*(-1+b)*p*p*pow(x10,b)*x20*(6*r*x10*x10+b*pow(x10,2*b)*x20*x20)))/(720*a*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20)+(1/(24*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20))*((x1-x10)*(x1-x10)*(x1-x10)*pow(x10,-3-4*b)*(4*b*(-1+p*p)*pow(x10,2*b)*x20*x20*(2*(1-2*b)*r*x10*x10+pow(x10,2*b)*x20*x20)-a*(-1+b)*p*pow(x10,1+b)*x20*(2*r*x10*x10+b*pow(x10,2*b)*x20*x20)-2*a*a*x10*x10*(2*r*x10*x10+b*(-1+p*p)*pow(x10,2*b)*x20*x20)))+((1/(24*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20))*(x1-x10)*(x1-x10)*(x2-x20)*(-2*a*a*a*(-1+p*p)*pow(x10,3+b)*x20-(-1+b)*b*p*(-1+p*p)*pow(x10,2*b)*x20*x20*(-6*r*x10*x10+(-4+b)*pow(x10,2*b)*x20*x20)+a*pow(x10,1+b)*x20*(-2*(8-5*p*p+b*(-14+11*p*p))*r*x10*x10+3*(-1+b)*b*p*p*pow(x10,2*b)*x20*x20)+2*a*a*p*x10*x10*(6*r*x10*x10+b*(-1+p*p)*pow(x10,2*b)*x20*x20)))/pow(x10,3*(1+b))+((1/(720*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20))*(x1-x10)*(x1-x10)*(x1-x10)*(x1-x10)*(7*a*a*a*a*x10*x10*x10*x10-30*b*(-1+p*p)*pow(x10,2*b)*x20*x20*(2*(1-4*b*b)*r*x10*x10+3*pow(x10,2*b)*x20*x20)+15*a*(-1+b)*b*p*pow(x10,1+b)*x20*(6*r*x10*x10+(2+b)*pow(x10,2*b)*x20*x20)+5*a*a*x10*x10*(12*(-1+4*b)*r*x10*x10+b*(4+7*b)*(-1+p*p)*pow(x10,2*b)*x20*x20)))/pow(x10,4*(1+b))+((1/(720*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20))*(x1-x10)*(x1-x10)*(x2-x20)*(x2-x20)*(2*a*a*a*(-23+44*p*p)*pow(x10,3+b)*x20-15*a*pow(x10,1+b)*x20*(2*(-16+13*p*p+b*(26-29*p*p))*r*x10*x10+3*(-1+b)*b*p*p*pow(x10,2*b)*x20*x20)-60*a*a*p*x10*x10*(12*r*x10*x10+b*(-1+p*p)*pow(x10,2*b)*x20*x20)-15*(-1+b)*b*p*pow(x10,2*b)*x20*x20*(6*(-3+2*p*p)*r*x10*x10+(-2+b)*(1+2*p*p)*pow(x10,2*b)*x20*x20)))/pow(x10,3*(1+b))+(1/(720*a*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20))*((x1-x10)*pow(x10,-3-2*b)*(x2-x20)*(x2-x20)*(x2-x20)*(240*a*a*(-2+5*p*p)*r*x10*x10*x10*x10-4*a*a*a*p*(-23+30*p*p)*pow(x10,3+b)*x20+30*(-1+b)*(-1+b)*b*p*p*pow(x10,4*b)*x20*x20*x20*x20+15*a*(-1+b)*p*pow(x10,1+b)*x20*(2*(19-12*p*p)*r*x10*x10+b*(1+2*p*p)*pow(x10,2*b)*x20*x20)))+((x1-x10)*(x2-x20)*(x2-x20)*(-4*a*(-5+8*p*p)*r*x10*x10*x10+4*a*a*p*(-1+p*p)*pow(x10,2+b)*x20-(-1+b)*p*pow(x10,b)*x20*(-6*(-3+2*p*p)*r*x10*x10+b*(1+2*p*p)*pow(x10,2*b)*x20*x20)))/(pow(x10,2*(1+b))*(24*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20))+((1/(720*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20*x20))*(x1-x10)*(x1-x10)*(x1-x10)*(x2-x20)*(-28*a*a*a*a*p*pow(x10,4+b)*x20+30*(-1+b)*b*p*(-1+p*p)*pow(x10,3*b)*x20*x20*x20*(-2*(1+b)*r*x10*x10+(-3+b)*pow(x10,2*b)*x20*x20)-30*a*b*pow(x10,1+2*b)*x20*x20*(-2*(6-3*p*p+b*(-10+7*p*p))*r*x10*x10+3*(-1+b)*p*p*pow(x10,2*b)*x20*x20)+60*a*a*a*x10*x10*x10*(4*r*x10*x10+b*(-1+p*p)*pow(x10,2*b)*x20*x20)-5*a*a*p*pow(x10,2+b)*x20*(18*(-1+5*b)*r*x10*x10+b*(-1+4*p*p+b*(-7+4*p*p))*pow(x10,2*b)*x20*x20)))/pow(x10,4*(1+b));

	double cOne=(12*r*r*x10*x10*x10*x10+12*(1-p*p+b*(-2+p*p))*r*pow(x10,2+2*b)*x20*x20+pow(x10,2*b)*x20*x20*(-4*a*a*(-1+p*p)*x10*x10+3*b*(2-b-2*p*p+2*b*p*p)*pow(x10,2*b)*x20*x20))/(pow(x10,2*(1+b))*(24*(-1+p*p)*x20*x20))+(pow(x10,-3-2*b)*(x2-x20)*(4*a*a*p*r*pow(x10,4+b)*x20-2*(-1+b)*(-1+b)*b*p*pow(x10,5*b)*x20*x20*x20*x20*x20-a*x10*(12*r*r*x10*x10*x10*x10-4*(-1+b)*r*pow(x10,2+2*b)*x20*x20-b*(6-3*b-4*p*p+4*b*p*p)*pow(x10,4*b)*x20*x20*x20*x20)))/(24*a*(-1+p*p)*x20*x20*x20)+((x1-x10)*pow(x10,-3-2*b)*(-4*a*a*r*x10*x10*x10*x10+2*a*(-1+b)*p*pow(x10,1+b)*x20*(-2*r*x10*x10+b*pow(x10,2*b)*x20*x20)+(-1+b)*(-12*r*r*x10*x10*x10*x10+b*(6-8*p*p+b*(-3+8*p*p))*pow(x10,4*b)*x20*x20*x20*x20)))/(24*(-1+p*p)*x20*x20)+((1/(1440*a*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20))*(x2-x20)*(x2-x20)*(-240*a*a*a*p*(-1+p*p)*r*pow(x10,5+b)*x20-4*a*a*a*a*(-1+p*p)*pow(x10,4+2*b)*x20*x20+15*(-1+b)*(-1+b)*p*p*pow(x10,2*b)*x20*x20*(2*r*x10*x10+b*pow(x10,2*b)*x20*x20)*(2*r*x10*x10+b*pow(x10,2*b)*x20*x20)-30*a*(-1+b)*p*pow(x10,1+b)*x20*(-12*r*r*x10*x10*x10*x10-4*b*r*pow(x10,2+2*b)*x20*x20+b*(6-6*p*p+b*(-5+6*p*p))*pow(x10,4*b)*x20*x20*x20*x20)+30*a*a*x10*x10*(4*(-5+8*p*p)*r*r*x10*x10*x10*x10-4*(-1+b)*(-1+p*p)*r*pow(x10,2+2*b)*x20*x20-b*(2-2*p*p*p*p+b*(-1+2*p*p*p*p))*pow(x10,4*b)*x20*x20*x20*x20)))/pow(x10,2*(2+b))+((1/(1440*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20))*(x1-x10)*(x1-x10)*(-4*a*a*a*a*(-1+p*p)*pow(x10,4+2*b)*x20*x20-30*a*(-1+b)*p*pow(x10,1+b)*x20*(-12*r*r*x10*x10*x10*x10-4*b*p*p*r*pow(x10,2+2*b)*x20*x20+b*(-2+b+2*p*p)*pow(x10,4*b)*x20*x20*x20*x20)+30*a*a*x10*x10*(12*r*r*x10*x10*x10*x10+4*(-1+2*b)*(-1+p*p)*r*pow(x10,2+2*b)*x20*x20-b*(2-2*p*p+b*(-1+2*p*p))*pow(x10,4*b)*x20*x20*x20*x20)+15*(-1+b)*pow(x10,2*b)*x20*x20*(4*(4-5*p*p+b*(-10+11*p*p))*r*r*x10*x10*x10*x10+4*(-1+b)*b*p*p*r*pow(x10,2+2*b)*x20*x20+b*(b*(-24+79*p*p-56*p*p*p*p)+12*(2-5*p*p+3*p*p*p*p)+b*b*(6-25*p*p+20*p*p*p*p))*pow(x10,4*b)*x20*x20*x20*x20)))/pow(x10,4*(1+b))+(1/(720*a*(-1+p*p)*(-1+p*p)*x20*x20*x20*x20))*(x1-x10)*pow(x10,-4-3*b)*(x2-x20)*(120*a*a*a*(-1+p*p)*r*pow(x10,5+b)*x20+4*a*a*a*a*p*(-1+p*p)*pow(x10,4+2*b)*x20*x20+30*a*(-1+b)*pow(x10,1+b)*x20*(4*(-5+2*p*p)*r*r*x10*x10*x10*x10-4*b*p*p*r*pow(x10,2+2*b)*x20*x20+b*(-3+4*p*p)*(2-2*p*p+b*(-1+2*p*p))*pow(x10,4*b)*x20*x20*x20*x20)+30*a*a*p*x10*x10*(-12*r*r*x10*x10*x10*x10+b*(3-3*p*p+b*(-2+3*p*p))*pow(x10,4*b)*x20*x20*x20*x20)-15*(-1+b)*(-1+b)*p*pow(x10,2*b)*x20*x20*(-6*(-1+b)*b*pow(x10,4*b)*x20*x20*x20*x20+p*p*(4*r*r*x10*x10*x10*x10+4*b*r*pow(x10,2+2*b)*x20*x20+b*(-6+7*b)*pow(x10,4*b)*x20*x20*x20*x20)));

	double cTwo=((1/(720*(-1+p*p)*x20*x20))*(4*a*a*a*a*(-1+p*p)*pow(x10,4+2*b)*x20*x20+15*(-1+b)*(-1+b)*pow(x10,2*b)*x20*x20*(-4*(-2+p*p)*r*r*x10*x10*x10*x10-4*b*p*p*r*pow(x10,2+2*b)*x20*x20+b*(12-6*b-12*p*p+11*b*p*p)*pow(x10,4*b)*x20*x20*x20*x20)+30*a*a*x10*x10*(4*r*r*x10*x10*x10*x10+b*(2-2*p*p+b*(-1+2*p*p))*pow(x10,4*b)*x20*x20*x20*x20)+120*a*(-1+b)*p*pow(x10,1+b)*x20*(2*r*r*x10*x10*x10*x10-b*r*pow(x10,2+2*b)*x20*x20+b*(2-2*p*p+b*(-1+2*p*p))*pow(x10,4*b)*x20*x20*x20*x20)))/pow(x10,2*(2+b));
	
	double cPoly=cMinusOne/Time + cZero + Time*cOne + 0.5*Time*Time*cTwo;
	double cDVP= sqrt(1-p*p)*a*pow(x1,b)*x2*x2*2*3.14159265358979323846*Time;

	return exp(cPoly)/cDVP;

}

double SABRProcess::InitsLow(double Time, double Deviations) const
{
	return LogAsset + (r - 0.5*Vol)*Time - Deviations*sqrt(Vol)*sqrt(Time);
}

double SABRProcess::InitsHigh(double Time, double Deviations) const
{
	return LogAsset + (r - 0.5*Vol)*Time + Deviations*sqrt(Vol)*sqrt(Time);
}


double SABRProcess::InitsLow(double sValue, double Time, double Deviations) const
{
	return sValue + (r - 0.5*Vol)*Time - Deviations*sqrt(Vol)*sqrt(Time);
}

double SABRProcess::InitsHigh(double sValue, double Time, double Deviations) const
{
	return sValue + (r - 0.5*Vol)*Time + Deviations*sqrt(Vol)*sqrt(Time);
}


double SABRProcess::InitvLow(double Time, double Deviations) const
{
	return InitvLow(Vol,Time,Deviations);
}

double SABRProcess::InitvHigh(double Time, double Deviations) const
{
	return InitvHigh(Vol,Time,Deviations);
}

double SABRProcess::InitvLow(double vValue, double Time, double Deviations) const
{
	double Guess = vValue - Deviations*sqrt(Time)*0.1;
	return max(Guess,0.001);
}

double SABRProcess::InitvHigh(double vValue, double Time, double Deviations) const
{
	double Guess = vValue + Deviations*sqrt(Time)*0.1;
	return Guess;
}


#endif
