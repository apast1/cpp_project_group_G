#ifndef GENERIC_QUAD_H
#define GENERIC_QUAD_H

#include "SVProcess.h"
#include "BetterGrid.h"
#include "SquareGrid2.h"
#include <omp.h>

const double DEVIATIONS = 5;
//heston
const double vDEVIATIONS = 2;
const double vsDEVIATIONS = 2;
//sabr
//const double vDEVIATIONS = 5;
const double TOLERANCE = 0.000000001;
const double TimeTOLERANCE = 0.1;

vector<double> SplitTime(double Time, unsigned N);

double SubIntegral(double sValue, double vValue, double Time, SVProcess& Proc, SquareGrid& GridPts);


double SubIntegral(double sValue, double vValue, double Time, SVProcess& Proc, SquareGrid& GridPts)
{
    //Get the parameters of the Grid
	double sLow = GridPts.GetsLow();
	double sHigh = GridPts.GetsHigh();
	double sStep = GridPts.GetsStep(); // NB: This is *half* of the "Simpson step"

	double vLow = GridPts.GetvLow();
	double vHigh = GridPts.GetvHigh();
	double vStep = GridPts.GetvStep();  // NB: This is *half* of the "Simpson step"

	// Compute the new integration range
	double NewsLow = max(Proc.InitsLow(sValue,Time,DEVIATIONS),sLow);
	double NewsHigh = min(Proc.InitsHigh(sValue,Time,DEVIATIONS),sHigh);
	double NewvLow = max(Proc.InitvLow(vValue,Time,vsDEVIATIONS),vLow);
	double NewvHigh = min(Proc.InitvHigh(vValue,Time,vsDEVIATIONS),vHigh);

	double Result=0.0;

	// If the integration range is empty, no point computing anything.
	if (NewsHigh > NewsLow && NewvHigh > NewvLow)
	{
		// Compute the entries of the grid corresponding to NewLow and NewHigh
		unsigned long isMinus = static_cast<unsigned long>(floor((NewsLow-sLow)/sStep));
		unsigned long isPlus = static_cast<unsigned long>(floor((NewsHigh-sLow)/sStep));

		unsigned long ivMinus = static_cast<unsigned long>(floor((NewvLow-vLow)/vStep));
		unsigned long ivPlus = static_cast<unsigned long>(floor((NewvHigh-vLow)/vStep));

		// Must check that the range has an odd number of values (for Simpson)

		if ((isPlus-isMinus)%2)
		{
			if (isPlus < GridPts.GetsSize()-1)
			{
				isPlus += 1;
			}
			else
			{
				isMinus -= 1;
			}
		}

		if ((ivPlus-ivMinus)%2)
		{
			if (ivPlus < GridPts.GetvSize()-1)
			{
				ivPlus += 1;
			}
			else
			{
				ivMinus -= 1;
			}
		}

        unsigned long iN = isPlus-isMinus+1;
        unsigned long vN = ivPlus-ivMinus+1;

        double Fact[iN][vN], Dent[iN][vN];

        # pragma omp parallel
        {
            # pragma omp for nowait
            for (unsigned long i = 0; i <iN ; i++)
            {
                for (unsigned long j = 0; j <vN; j++)
                {
                    if ((i == 0 || i == iN-1)&&(j == 0 || j == vN-1))
                    {
                        Fact[i][j]=1.0;

                    }
                    else if (i == 0 || i == iN-1)
                    {
                        if(j%2)
                        {
                            Fact[i][j]=4.0;
                        }
                        else
                        {
                            Fact[i][j]=2.0;
                        }
                    }
                    else if (j == 0 || j == vN-1)
                    {
                        if(i%2)
                        {
                        Fact[i][j]=4.0;
                        }
                        else
                        {
                        Fact[i][j]=2.0;
                        }
                    }
                    else
                    {
                        if ((i%2)&&(j%2))
                        {
                            Fact[i][j]=16.0;
                        }
                        else if (!(i%2)&&!(j%2))
                        {
                            Fact[i][j]=4.0;
                        }
                        else
                        {
                            Fact[i][j]=8.0;
                        }
                    }
                }
            }

            # pragma omp for collapse(2)
            for (unsigned long i = 0; i <iN ; i++)
            {
                for (unsigned long j = 0; j <vN; j++)
                {
                    Dent[i][j]=GridPts.GetValue(i+isMinus,j+ivMinus)*Proc.Density(Time,sLow+(i+isMinus)*sStep,vLow+(j+ivMinus)*vStep,sValue,vValue);
                }
            }

            # pragma omp for reduction(+:Result) nowait collapse(2)
            for (unsigned long i = 0; i <iN ; i++)
            {
                for (unsigned long j = 0; j <vN; j++)
                {
                    Result += Fact[i][j]*Dent[i][j];
                }
            }
        }
        Result *= exp(-Time*Proc.GetRate())*vStep*sStep/9.0;
	}
    return Result;
}

vector<double> SplitTime(double Time, unsigned N)
{
	double tmp = Time/static_cast<double>(N);

	vector<double> Result;

	for (unsigned i = 1; i < N+1; i++)
	{
		Result.push_back(i*tmp);
	}

	return Result;
};

#endif
