#ifndef MULTILEVEL_H
#define MULTILEVEL_H

#include <iostream>

#include "Processes.h"
#include "MersenneTwister.h"
#include "Option.h" 

using std::cout; 

class MLData
{
public:
	MLData(double Eps_, unsigned long Divisor_);
	bool DoneYet() const;
	double Bound() const;
	double GetSum() const;
	double GetVariance(unsigned long l) const; 
	void UpdateEvolutions(Process& Proc, Option& Opt); 
	void UpdateNewNs(double Time); 
	void UpdateFirstEvo(Process& Proc, Option& Opt);
	void UpdateLaterEvo(Process& Proc, Option& Opt, unsigned long L); 
	void NewN(unsigned long N); 
private:
	double Eps;
	unsigned long Divisor; 
	vector<double> Sums;
	vector<double> Squares;
	vector<unsigned long> OldNs;
	vector<unsigned long> NewNs;
	MTRand MT; 
};

double MLPriceOption(Process& Proc, Option& Opt, double Eps, unsigned long Divisor);

vector<double> GetTimeSteps(vector<double> Times); 


////////////////////////////////////////////////////////////////////////////////////////////////


double MLPriceOption(Process& Proc, Option& Opt, double Eps, unsigned long Divisor)
{
	MLData ML(Eps,Divisor); 

	unsigned long L = 0; 

	do
	{	
		if (L == 0)
		{
			ML.NewN(10000);
		}
		else
		{
			ML.NewN(2500);
		}
		ML.UpdateEvolutions(Proc,Opt);
		ML.UpdateNewNs(Opt.GetExpiry());
		ML.UpdateEvolutions(Proc,Opt); 
		L++; 
	}
	while (ML.DoneYet() == 0);

	return exp(-Proc.GetRate()*Opt.GetExpiry())*ML.GetSum(); 
}


MLData::MLData(double Eps_, unsigned long Divisor_) : Eps(Eps_), Divisor(Divisor_)
{
	vector<double> Sums;
	vector<double> Squares;
	vector<unsigned long> OldNs;
	vector<unsigned long> NewNs; 
	MTRand MT; 
}
 
bool MLData::DoneYet() const
// Checks whether the multilevel result is accurate yet
{
	if (Sums.size() < 2)
	{
		return 0;
	}
	else
	{
		unsigned long L = Sums.size(); 
		double Ya = fabs(Sums[L-2])/static_cast<double>(NewNs[L-2]);
		double Yb = fabs(Sums[L-1])/static_cast<double>(NewNs[L-1]); 
		double m = max(Ya/static_cast<double>(Divisor),Yb); 
		if (m < Bound())
		{
			return 1;
		}
		else
		{
			return 0;
		}
		return 1; 
	}
}

double MLData::Bound() const
// The bound for the convergence check of multilevel routines
{
	return (Divisor-1.0)*Eps/sqrt(2.0); 
}

double MLData::GetSum() const
// Sums over the "Sums", i.e. returns the estimate for E[PayOff]
{
	double Result = 0;

	for (unsigned long i = 0; i < Sums.size(); i++)
	{
		Result += Sums[i]/static_cast<double>(OldNs[i]);
	}

	return Result;
}

void MLData::UpdateNewNs(double Time)
{
	double Konst = 0;
	unsigned long L = OldNs.size(); 

	for (unsigned long i = 0; i < L; i++)
	{
		Konst += sqrt(GetVariance(i)*pow(static_cast<double>(Divisor),static_cast<int>(i))); 
	}

	Konst *= 2.0/(Eps*Eps); 

	for (unsigned long i = 0; i < L; i++)
	{
		double Ndouble = ceil(Konst*sqrt(GetVariance(i)/pow(static_cast<double>(Divisor),static_cast<int>(i))));
		NewNs[i] = static_cast<unsigned long>(Ndouble); 
	}
}

void MLData::NewN(unsigned long N) 
{
	NewNs.push_back(N);
	OldNs.push_back(0); 
}

void MLData::UpdateEvolutions(Process& Proc, Option& Opt)
{
	unsigned long L = NewNs.size(); 

	for (unsigned long k = 0; k < L; k++)
	{
		if (k == 0)
		{
			UpdateFirstEvo(Proc,Opt); 
		}
		else
		{
			UpdateLaterEvo(Proc,Opt,k);
		}
	}
}

void MLData::UpdateFirstEvo(Process& Proc, Option& Opt)
{
	if (Sums.empty() || Squares.empty())
	{
		Sums.push_back(0);
		Squares.push_back(0);
	}

	if (NewNs[0] > OldNs[0])
	{
		unsigned long Paths = NewNs[0] - OldNs[0];
		vector<double> TimeSteps = Opt.GetTimeSteps(); 
		for (unsigned i = 0 ; i < Paths; i++)
		{
				Process* NewProc = Proc.clone();

				// Simulate a sample path
				vector<double> Obs;
				for (unsigned long k = 0; k < TimeSteps.size(); k++)
				{
					// Evolve the process over the kth timestep and store the value in Obs
					NewProc->Evolve(TimeSteps[k],MT.randNorm()); 
					Obs.push_back(NewProc->GetValue()); 
				}
				//Then Compute the payoff
				double p = Opt.PayOff(Obs); 
	
				Sums[0] += p; 
				Squares[0] += p*p; 
				delete NewProc;
		}
		OldNs[0] = NewNs[0]; 
	}
}

void MLData::UpdateLaterEvo(Process& Proc, Option& Opt, unsigned long l)
{ 
	if (NewNs[l] > OldNs[l])
	{
		unsigned long Paths = NewNs[l] - OldNs[l]; 
		vector<double> TimeSteps = Opt.GetTimeSteps(); 
		
		if (Sums.size() < l+1)
		{
			Sums.push_back(0);
			Squares.push_back(0); 
		}

		for (unsigned long PathNumber = 0; PathNumber < Paths; PathNumber++)
		{
			vector<double> FineObservations;
			vector<double> CoarseObservations; 

			Process* FineProc = Proc.clone(); 
			Process* CoarseProc = Proc.clone(); 
	
			for (unsigned long StepNumber =0 ; StepNumber < TimeSteps.size(); StepNumber++)
			{
				double SmallTimeStep = TimeSteps[StepNumber]/pow(static_cast<double>(Divisor),static_cast<int>(l)); 
				double BigTimeStep = static_cast<double>(Divisor)*SmallTimeStep; 
				unsigned long NumberOfSteps = Power(Divisor,l); 

				for (unsigned k = 0; k < NumberOfSteps; k++)
				// Each step of the path simulation done in NumberOfSteps parts
				{
					double CoarseRand = 0; 
	
					for (unsigned j = 0; j < Divisor; j++)
					// Evolve Proc in M steps & compute the sum of the random draws
					{
						double r = MT.randNorm();
						CoarseRand += r; 
						FineProc->Evolve(SmallTimeStep,r);
					}
					CoarseRand = CoarseRand/sqrt(static_cast<double>(Divisor));
					CoarseProc->Evolve(BigTimeStep,CoarseRand); 
				}

				FineObservations.push_back(FineProc->GetValue()); 
				CoarseObservations.push_back(CoarseProc->GetValue()); 
			}
			double x = Opt.PayOff(FineObservations) - Opt.PayOff(CoarseObservations);
			Sums[l] += x;
			Squares[l] += x*x; 
			delete FineProc;
			delete CoarseProc;
		}

	OldNs[l] += Paths; 
	}
}

double MLData::GetVariance(unsigned long l) const
{
	double Nl = static_cast<double>(OldNs[l]);
	return (Squares[l]/Nl - Sums[l]*Sums[l]/(Nl*Nl)); 
}


#endif