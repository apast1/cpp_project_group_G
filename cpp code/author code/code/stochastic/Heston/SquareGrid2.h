#ifndef SQUARE_GRID_H
#define SQUARE_GRID_H

#include <fstream>

#include "BetterGrid.h"

class SquareGrid
{
public:	
	SquareGrid();
	SquareGrid(const SquareGrid& Orig); 
	SquareGrid(double sLow_, double vLow_, double sStep_, double vStep_, unsigned long sN_, unsigned long vN_); 
	SquareGrid(double sLow_, double vLow_, double sStep_, double vStep_, unsigned long sN_, unsigned long vN_, double **Values_); 
	SquareGrid& operator=(const SquareGrid& Old); 
	SquareGrid* clone() const;
	double GetsLow() const;
	double GetvLow() const;
	double GetsStep() const;
	double GetvStep() const;
	double GetsHigh() const;
	double GetvHigh() const;
	unsigned long GetsSize() const;
	unsigned long GetvSize() const;
	double GetValue(int x, int y) const;
	void SetValue(unsigned long x, unsigned long y, double New); 
	void PrintGrid() const;
	void PrintIntoFile(char* Name);
	~SquareGrid(); 
private:
	double sLow;
	double vLow;
	double sStep;
	double vStep;
	unsigned long sN;
	unsigned long vN;
	double **Values;
};


SquareGrid::SquareGrid()
{
	sLow = 0.0;
	vLow = 0.0;
	sStep = 0.0;
	vStep = 0.0;
	sN = 0;
	vN = 0;
	Values = NULL;
}

SquareGrid::SquareGrid(double sLow_, double vLow_, double sStep_, double vStep_, unsigned long sN_, unsigned long vN_) :
sLow(sLow_), vLow(vLow_), sStep(sStep_), vStep(vStep_), sN(sN_), vN(vN_)
{
	Values = new double*[sN+1];
	for (unsigned long i = 0; i < sN+1; i++)
	{
		Values[i] = new double[vN+1];
	}
}

SquareGrid::SquareGrid(double sLow_, double vLow_, double sStep_, double vStep_, unsigned long sN_, unsigned long vN_, double **Values_) : 
sLow(sLow_), vLow(vLow_), sStep(sStep_), vStep(vStep_), sN(sN_), vN(vN_) 
{
	Values = new double*[sN+1];
	for (unsigned long i = 0; i < sN+1; i++)
	{
		Values[i] = new double[vN+1];
		memcpy(Values[i],Values_[i],(vN+1)*sizeof(double)); 
	}
}

SquareGrid::SquareGrid(const SquareGrid& Orig)
{
	sLow = Orig.sLow;
	vLow = Orig.vLow;
	sStep = Orig.sStep;
	vStep = Orig.vStep;
	sN = Orig.sN;
	vN = Orig.vN;

	Values = new double*[sN+1];
	for (unsigned long i = 0; i < sN+1; i++)
	{
		Values[i] = new double[vN+1];
		memcpy(Values[i],Orig.Values[i],(vN+1)*sizeof(double));
	}
}

SquareGrid& SquareGrid::operator=(const SquareGrid& Orig)
{
	if (this != &Orig)
	{
		for (unsigned long i = 0; i < sN+1; i++)
		{
			delete [] Values[i];
		}
		delete [] Values;

		sLow = Orig.sLow;
		vLow = Orig.vLow;
		sN = Orig.sN;
		sStep = Orig.sStep;
		vN = Orig.vN;
		vStep = Orig.vStep;

		Values = new double*[sN+1];

		for (unsigned long i = 0; i < sN+1; i++)
		{
			Values[i] = new double[vN+1];
			memcpy(Values[i],Orig.Values[i],(vN+1)*sizeof(double));
		}
	}

	return *this; 
}

double SquareGrid::GetsLow() const
{
	return sLow;
}

double SquareGrid::GetsHigh() const
{
	return sLow + sN*sStep;
}

double SquareGrid::GetsStep() const
{
	return sStep;
}


unsigned long SquareGrid::GetsSize() const
{
	return sN+1;
}

double SquareGrid::GetvLow() const
{
	return vLow;
}

double SquareGrid::GetvHigh() const
{
	return vLow + vN*vStep;
}

double SquareGrid::GetvStep() const
{
	return vStep;
}

unsigned long SquareGrid::GetvSize() const
{
	return vN+1;
}

double SquareGrid::GetValue(int x, int y) const
{
	return Values[x][y];
}

void SquareGrid::SetValue(unsigned long x, unsigned long y, double New) 
{
	Values[x][y] = New;
}


void SquareGrid::PrintGrid() const
{
	for (unsigned long i = 0; i < sN+1; i++)
	{
		for (unsigned long j = 0; j < vN+1; j++)
		{
			std::cout << Values[i][j] << " ";
		}
		std::cout << "\n"; 
	}
}

SquareGrid* SquareGrid::clone() const
{
	return new SquareGrid(*this);
}

SquareGrid::~SquareGrid()
{
	for (unsigned long i = 0; i < sN+1; i++)
	{
		delete [] Values[i];
	}
	delete [] Values;
}

void SquareGrid::PrintIntoFile(char* Name)
{
	std::ofstream filu;
	filu.open(Name);

	for (unsigned long i = 0; i < sN+1; i++)
	{
		for (unsigned long j = 0; j < vN+1; j++)
		{
			filu << Values[i][j] << ",";
		}
		filu << "\n";
	}

	filu.close();
}

#endif
