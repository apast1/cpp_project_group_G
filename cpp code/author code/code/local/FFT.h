#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex> 
#include <fstream>

#include "Processes.h"
#include "MiscMaths.h"
#include "Constants.h"

using std::complex; 
using std::vector; 
using std::fstream; 
using namespace std; 

typedef complex<double> Complex; 


/////////////////////////////////////////////////////////////////////////////////////////////////////

class FTData
{
	// Stores the data needed for working with characteristic function chi and its Fourier transform f
	// Write s for the variable in the transform and xi for the variable of the characteristic function
	// N = (# of the values of chi(xi) to be given) = (# of the values f(s) to be obtained) 
	// sLow = Lowest value of s to be obtained
	// Bound = The upper bound of the Fourier integral "\int_0^B exp(...) chi(xi) dxi"
	// GetsStep and GetxiStep return 2*pi/Bound and Bound/N, respectively.
public:
	FTData(unsigned long N_, double sLow_, double Bound_);
	unsigned long GetN() const; 
	double GetsLow() const; 
	double GetBound() const; 
	double GetsStep() const; 
	double GetxiStep() const;
private:
	unsigned long N; 
	double sLow;
	double Bound;
};

vector<Complex> FFT(vector<Complex> Vec);
vector<Complex> InvFFT(vector<Complex> Vec); 

vector<Complex> MakeCharVector(double Time, Process& Proc, FTData& Data);
vector<double> GetDensityFromChar(double Time, Process& Proc, FTData& Data); 

//---------------------------------------------------------------------------------------------------

void WriteIntoFile(vector<double> Vec, Process& Proc, double Time, double B); 

double CorrW(double theta);

Complex CorrA(double theta);

vector<Complex> TakeEvens(vector<Complex> Vec); 
vector<Complex> TakeOdds(vector<Complex> Vec); 

/////////////////////////////////////////////////////////////////////////////////////////////////////

FTData::FTData(unsigned long N_, double sLow_, double Bound_) :
	N(N_), Bound(Bound_), sLow(sLow_)
{
}

unsigned long FTData::GetN() const
{
	return N;
}

double FTData::GetBound() const
{
	return Bound;
}

double FTData::GetsLow() const
{
	return sLow ;
}

double FTData::GetxiStep() const
// Given the integration Bound and the number of Steps, returns the step size xi
{
	return Bound/static_cast<double>(N);
}

double FTData::GetsStep() const
{
	return 2.0*PI/Bound; 
}
vector<Complex> FFT(vector<Complex> Vec)
// Computer the Fast Fourier Transform of a Complex vector using Cooley-Tukey algorithm
// We normalise discrete Fourier transform so that the exponential term is exp(-2pi*j*k/N)
{
	vector<Complex> Result;

	if (Vec.size() == 2)
	{
		// If size = 2, compute by hand
		Result.push_back(Vec[0]+Vec[1]);
		Result.push_back(Vec[0]-Vec[1]);
	}
	else
	{
		// If size > 2, we split the vector into two vectors of length size/2,
		// according to the Cooley-Tukey algo

		unsigned long M = Vec.size()/2; 

		// Precompute 2pi*i/N; 
		Complex Argu = -(2.0/static_cast<double>(Vec.size()))*PI*I;

		// Form the shorter vectors consisting of Even/Odd entries in Vec
		// and Fourier transform them
		vector<Complex> EvensFT=FFT(TakeEvens(Vec));
		vector<Complex> OddsFT=FFT(TakeOdds(Vec)); 

		// Compute the DFT of Vec by putting together the shorted transforms
		for (unsigned long k = 0; k < M; k++)
		{
			Complex Value = EvensFT[k] + exp(static_cast<double>(k)*Argu)*OddsFT[k]; 
			Result.push_back(Value); 
		}
		for (unsigned long k = M ; k < 2*M ; k++)
		{
			Complex Value = EvensFT[k-M] - exp(static_cast<double>(k-M)*Argu)*OddsFT[k-M]; 
			Result.push_back(Value);
		}
	}

	return Result; 
}

vector<Complex> InvFFT(vector<Complex> Vec)
{
	vector<Complex> Result;
	Complex I(0.0,1.0); 
	if (Vec.size() == 2)
	{
		Result.push_back(0.5*(Vec[0]+Vec[1]));
		Result.push_back(0.5*(Vec[0]-Vec[1]));
	}
	else
	{
		unsigned long M = Vec.size()/2; 
	
		Complex Argu = (2/static_cast<double>(Vec.size()))*PI*I;

		vector<Complex> EvensFT=InvFFT(TakeEvens(Vec));
		vector<Complex> OddsFT=InvFFT(TakeOdds(Vec)); 

		for (unsigned long i = 0; i < M; i++)
		{
			Complex Value = EvensFT[i] + exp(-static_cast<double>(i)*Argu)*OddsFT[i]; 
			Result.push_back(0.5*Value); 
		}
		for (unsigned long i = M ; i < 2*M ; i++)
		{
			Complex Value = EvensFT[i-M] - exp(-static_cast<double>(i-M)*Argu)*OddsFT[i-M]; 
			Result.push_back(0.5*Value);
		}
	}

	return Result; 
}

vector<Complex> MakeCharVector(double Time, Process& Proc, FTData& Data)
// Prepares a vector for Fourier transform by computing the appropriate values of 
// the characteristic function and by giving these values weights according to 
// Simpson's Rule

// The point is to integrate exp(i*xi*s) chi(xi) from 0 to a Bound given in the FTData
// This does _not_ attempt to integrate from -infty to infty
{
	vector<Complex> Result;

	unsigned long N = Data.GetN(); 
	double xiStep = Data.GetxiStep(); 
	double sLow = Data.GetsLow(); 

	double xi = 0; 

	// Enter the value of chi at 0 together with factor from Simpson's rule
	Complex NewEntry = Proc.Char(Time,xi)*xiStep/3.0; 
	Result.push_back(NewEntry);

	// Then the rest
	for (unsigned long i = 1; i < N; i++)
	{
		xi = static_cast<double>(i)*xiStep; 

		// First a fudge term needed to turn the integral into FFT'able data
		// times a factor from Simpson's
		NewEntry = exp(-I*sLow*xi)*xiStep/3.0; 

		// Then the values chi(xi)
		NewEntry *= Proc.Char(Time,xi); 

		// Finally the weights according to Simpson's rule
		if (i%2 == 0)
		{
			NewEntry *= 2.0;
		}
		else
		{
			NewEntry *= 4.0; 
		}

		Result.push_back(NewEntry);
	}

	return Result;
}

vector<double> GetDensityFromChar(double Time, Process& Proc, FTData& Data)
// Computes the density function of a process as the DFT of the characteristic function
{
	// Prepare the data for FFT
	vector<Complex> CharVector = MakeCharVector(Time,Proc,Data); 

	// Do the FFT and take the real part 
	vector<Complex> FT = FFT(CharVector); 
	vector<double> Result = realpart(FT);

	// Divide by pi
	Multiply(Result,OneByPI); 

	return Result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

// Some auxiliary algorithms

vector<Complex> TakeEvens(vector<Complex> Vec)
{
	vector<Complex> Result;
	unsigned long i = 0; 
	while (i < Vec.size())
	{
		Result.push_back(Vec[i]);
		i = i+2; 
	}

	return Result;
}

void WriteIntoFile(vector<double> Vec, Process& Proc, double Time, double B)
{
	fstream filu("Numeroita.hjh",ios::out);

	unsigned long N = Vec.size();

	double omega=0; 
	unsigned long i = 0; 

	do
	{
		filu << omega << " " << Vec[i] << " " << Proc.Density(Time,omega) << "\n"; 
		omega += PI/(2.0*B); 
		i++; 
	}
	while (omega < 10); 

	filu.close(); 
}


vector<Complex> TakeOdds(vector<Complex> Vec)
{
	vector<Complex> Result;
	unsigned long i = 1;
	do 
	{
		Result.push_back(Vec[i]);
		i = i+2;
	}
	while (i < Vec.size()); 
	return Result;
}

double CorrW(double theta)
{
	return 2*(1-cos(theta))/(theta*theta); 
}

Complex CorrA(double theta)
{
	Complex Result(-(1-cos(theta))/(theta*theta), (theta-sin(theta))/(theta*theta));
	return Result;
}



#endif