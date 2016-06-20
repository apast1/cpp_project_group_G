#ifndef GRID_H
#define GRID_H

// const unsigned long Steps = 10; 
class Grid
{
public:
	Grid();
	Grid(const Grid& Orig); 
	Grid(double Low_, double Step_, unsigned long Steps_); 
	Grid(double Low_, double Step_, unsigned long Steps_, double *Values_); 
	Grid& operator=(const Grid& Old); 
	Grid* clone() const;
	double GetLow() const;
	double GetStep() const;
	double GetNthValue(int N) const;
	double GetHigh() const;
	double* GetVals() const;
	unsigned long GetSize() const;
	void SetLow(double NewLow);
	void SetStep(double NewStep);
	void SetNthValue(unsigned long N, double New); 
	~Grid(); 
	void PrintGrid() const;
	void PrintSteps() const;
private:
	double Low; 
	double Step;
	unsigned long Steps;
	double *Values;
};

Grid::Grid() : Low(0),Step(0),Steps(1)
{
	Values = NULL;
}

Grid::Grid(const Grid& Orig)
{
	Low = Orig.Low;
	Step = Orig.Step;
	Steps = Orig.Steps;
	Values = new double[Steps];
	memcpy(Values,Orig.Values,Steps*sizeof(double)); 
}

Grid::Grid(double Low_, double Step_, unsigned long Steps_, double *Values_) :
Low(Low_), Step(Step_), Steps(Steps_)
{
	Values = new double[Steps]; 
	memcpy(Values,Values_,Steps_*sizeof(double)); 
}

Grid::Grid(double Low_, double Step_,unsigned long Steps_) :
Low(Low_), Step(Step_),Steps(Steps_)
{
	Values = new double[Steps_];
}

Grid& Grid::operator=(const Grid& original)
{
	if (this != &original)
	{
		Low = original.Low;
		Step = original.Step;
		Steps = original.Steps;
		delete [] Values;
		Values = new double[Steps];
		memcpy(Values,original.Values,Steps*sizeof(double)); 
	}
	return *this;
}

Grid* Grid::clone() const
{
	return new Grid(*this);
}

double* Grid::GetVals() const
{
	return Values;
}

unsigned long Grid::GetSize() const
{
	return Steps;
}

void Grid::SetStep(double NewStep)
{
	Step = NewStep;
}

void Grid::SetLow(double NewLow)
{
	Low = NewLow;
}

void Grid::SetNthValue(unsigned long N, double New)
{
	Values[N] = New;
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
	return static_cast<double>(Low+(Steps-1)*Step); 
}

Grid::~Grid()
{
		delete [] Values;
}

void Grid::PrintGrid() const
{
	for (unsigned long i = 0; i < Steps; i++)
	{
		std::cout << Values[i] << " "; 
	}
}

void Grid::PrintSteps() const
{
	for (unsigned long i = 0; i < Steps; i++)
	{
		std::cout << Low+i*Step << " ";
	}
}

#endif