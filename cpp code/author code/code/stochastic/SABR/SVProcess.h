#ifndef SVPROCESS_H
#define SVPROCESS_H


class SVProcess
{
public:

	virtual double InitsLow(double Time, double Deviations) const = 0;
	virtual double InitsHigh(double Time, double Deviations) const = 0;
	virtual double InitvLow(double Time, double Deviations) const = 0;
	virtual double InitvHigh(double Time, double Deviations) const = 0;
	virtual double InitsLow(double sValue, double Time, double Deviations) const=0;
	virtual double InitsHigh(double sValue, double Time, double Deviations) const=0;
//	virtual double InitsLow(double sValue, double vValue, double Time, double Deviations) const=0;
//	virtual double InitsHigh(double sValue, double vValue, double Time, double Deviations) const=0;
	virtual double InitvLow(double vValue, double Time, double Deviations) const=0;
	virtual double InitvHigh(double vValue, double Time, double Deviations) const=0;
	virtual double Density(double Time, double AssetValue, double VarValue) const = 0;
	virtual double Density(double Time, double AssetValue, double VarValue, double AssetBase, double VarBase) const = 0;
	virtual double GetValue() const = 0;
	virtual double GetRate() const = 0;
	virtual double GetVar() const = 0;
	virtual void SetValues(double NewValue, double NewVar) = 0; 
}; 

#endif
