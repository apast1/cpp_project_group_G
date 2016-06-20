#ifndef OPTION_H
#define OPTION_H


class Option
{
public:
	virtual double PayOff(vector<double> Observations) const=0;
	virtual double GetExpiry() const=0; 
	virtual vector<double> GetTimeSteps() const=0;
}; 

class EuroCall : public Option
{
public:
	EuroCall(double Strike_, double Time_); 
	double PayOff(vector<double> Observations) const;
	double PayOff(double Value) const;
	double GetExpiry() const;
	vector<double> GetTimeSteps() const;
private:
	double Strike;
	double Time; 
};

class AsianCall : public Option
{
public:
	AsianCall(vector<double> Strikes, vector<double> Times); 
	double PayOff(vector<double> Observations) const;
	double GetExpiry() const;
	vector<double> GetTimeSteps() const;
private:
	vector<double> Strikes;
	vector<double> Times;
}; 

class LookBackCall : public Option
{
public:
	LookBackCall(vector<double> Strikes, vector<double> Times); 
	double PayOff(vector<double> Observation) const;
	double GetExpiry() const;
	vector<double> GetTimeSteps() const;
private:
	vector<double> Strikes;
	vector<double> Times;
}; 

vector<double> ComputeSteps(vector<double>);

///////////////////////////////////////////////////////////////////////////////////////

EuroCall::EuroCall(double Strike_, double Time_) : Strike(Strike_), Time(Time_)
{
}

double EuroCall::PayOff(double Value) const
{
	return max(Value-Strike,0.0);
}

double EuroCall::PayOff(vector<double> Values) const
{
	return max(Values[0]-Strike,0.0);
}

double EuroCall::GetExpiry() const
{
	return Time;
}

vector<double> EuroCall::GetTimeSteps() const
{
	vector<double> Result;
	Result.push_back(Time);
	return Result;
}

AsianCall::AsianCall(vector<double> Strikes_, vector<double> Times_) 
: Strikes(Strikes_), Times(Times_)
{
}

double AsianCall::PayOff(vector<double> Observations) const
{
	double Result = 0; 
	for (unsigned long i = 0; i < Observations.size(); i++)
	{
		Result += max(Observations[i]-Strikes[i],0.0); 
	}

	return Result/static_cast<double>(Observations.size());
}

double AsianCall::GetExpiry() const
{
	return Times.back(); 
}

vector<double> AsianCall::GetTimeSteps() const
{
	return ComputeSteps(Times); 
}

LookBackCall::LookBackCall(vector<double> Strikes_, vector<double> Times_) 
: Strikes(Strikes_), Times(Times_)
{
}

double LookBackCall::PayOff(vector<double> Observations) const
{
	double Result = 0;

	for (unsigned long i = 0; i < Observations.size(); i++)
	{
		double x = max(Observations[i]-Strikes[i],0.0);
		if (x > Result)
		{
			Result = x;
		}
	}

	return Result; 
}

double LookBackCall::GetExpiry() const
{
	return Times.back();
}

vector<double> LookBackCall::GetTimeSteps() const
{
	return ComputeSteps(Times);
}

vector<double> ComputeSteps(vector<double> Times)
{
	vector<double> Result;
	Result.push_back(Times[0]); 

	for (unsigned long i = 0; i < Times.size() - 1; i++)
	{
		Result.push_back(Times[i+1]-Times[i]); 
	}

	return Result; 
}


#endif