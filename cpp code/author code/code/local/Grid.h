#ifndef GRID_H
#define GRID_H

#include <vector>

// Stores values of functions, as needed for QUAD routines
// The idea is that Values[i] = f(Low+i*Step)

using std::vector; 

class Grid
{
public: 
	Grid(double Low_, double Step_, vector<double> Values_); 
	double GetLow() const;
	double GetHigh() const;
	double GetStep() const;
	double GetNthValue(int N) const;
	unsigned long GetSize() const;
private:
	double Low;
	double Step;
	vector<double> Values;
};

Grid::Grid(double Low_, double Step_, vector<double> Values_) :
Low(Low_), Step(Step_), Values(Values_)
{
}

unsigned long Grid::GetSize() const
{
	return Values.size();
}

double Grid::GetLow() const
{
	return Low;
}

double Grid::GetStep() const
{
	return Step;
}

double Grid::GetNthValue(int N) const
{
	return Values[N];
}

double Grid::GetHigh() const
{
	return static_cast<double>(Low+(Values.size()-1)*Step); 
}

#endif