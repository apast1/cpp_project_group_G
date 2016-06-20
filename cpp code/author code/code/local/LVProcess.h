#ifndef LVPROCESS_H
#define LVPROCESS_H

class LVProcess
{
public:
	virtual void Evolve(double Time, MTRand& MT) = 0; 
	virtual void Evolve(double Time, unsigned long N, MTRand& MT) = 0;
	// double InitVol() const;
	// double EuroPut(double Time, double Strike, unsigned long Steps) const; 
	virtual double InitLow(double Time, double Deviations) const = 0;
	virtual double InitHigh(double Time, double Deviations) const = 0;
	virtual double InitLow(double Base, double Time, double Deviations) const = 0;
	virtual double InitHigh(double Base, double Time, double Deviations) const = 0;
	// double EuroPutViaAppr(double Time, double Strike, unsigned long Steps) const; 
	// double EuroCallViaAppr(double Time, double Strike, unsigned long Steps) const;
	virtual double Density(double Time, double Value) const = 0;
	virtual double Density(double Time, double Value, double Base) const = 0;
	// double ApprZeroDensity(double Time, double Value, double Base) const;
	// double ApprZeroDensity(double Time, double Value) const;
	// double ApprDensity(double Time, double Value, double Base) const;
	// double ApprDensity(double Time, double Value) const;
	// double ZeroProb(double Time, unsigned long N) const;
	virtual double GetValue() const = 0;
	virtual double GetRate() const = 0;
	virtual void SetValue(double NewValue) = 0; 
}; 

#endif