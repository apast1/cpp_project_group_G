#include <cmath>
#include <complex>

#include "Constants.h" 

using std::complex; 

typedef complex<double> Complex;

class CxVec
{
public:
	CxVec(const CxVec& Orig); 
	CxVec(unsigned long Size_);
	CxVec(Complex *Vec_, unsigned long Size_);
	~CxVec();
	CxVec clone() const;
	CxVec fft() const;
	CxVec TakeOdds() const;
	CxVec TakeEvens() const;
	CxVec& operator=(const CxVec& Other); 
	Complex GetEntry(unsigned long Index) const;
	void SetEntry(unsigned long Index, Complex Value);
	unsigned long GetSize() const; 
private:
	unsigned long Size; 
	Complex *Vec;
};


CxVec::CxVec(const CxVec& Orig)
{
	Size = Orig.Size;
	Vec = new Complex[Size];
	memcpy(Vec,Orig.Vec,Size*sizeof(Complex)); 
}

CxVec::CxVec(unsigned long Size_) : Size(Size_) 
{
	Vec = new Complex[Size];
}



CxVec::~CxVec()
{
	delete [] Vec;
}


CxVec::CxVec(Complex *Vec_, unsigned long Size_) : Vec(Vec_), Size(Size_)
{
}

CxVec& CxVec::operator=(const CxVec& Other)
{
	if (this != &Other)
	{
		Size = Other.GetSize();
		Complex *NewArray = new Complex[Size]; 
		delete [] Vec; 

		for (unsigned long i = 0; i < Size; i++)
		{
			NewArray[i] = Other.GetEntry(i); 
		}

		Vec = NewArray; 
	}
	return *this; 
}

/*
CxVec CxVec::clone() const
{
	Complex *NewVec = new Complex[Size];
	for (unsigned long i = 0; i < Size; i++)
	{
		NewVec[i] = Vec[i];
	}
	CxVec Result(NewVec,Size);
	return Result;
}
*/

void CxVec::SetEntry(unsigned long Index, Complex Value)
{
	Vec[Index] = Value;
}

unsigned long CxVec::GetSize() const
{
	return Size; 
}

Complex CxVec::GetEntry(unsigned long Index) const
{
	return Vec[Index];
}


CxVec CxVec::TakeEvens() const
{
	Complex *NewVec = new Complex[Size/2]; 

	for (unsigned long i = 0; i < Size/2; i++)
	{
		NewVec[i] = Vec[2*i];
	}

	CxVec Result(NewVec,Size/2);
	return Result;
}

CxVec CxVec::TakeOdds() const
{
	Complex *NewVec = new Complex[Size/2]; 

	for (unsigned long i = 0; i < Size/2; i++)
	{
		NewVec[i] = Vec[2*i+1];
	}

	CxVec Result(NewVec,Size/2);
	return Result;
}

CxVec CxVec::fft() const
{
	if (Size == 2)
	{
		Complex *NewVec = new Complex[2];
		NewVec[0] = Vec[0] + Vec[1];
		NewVec[1] = Vec[0] - Vec[1];
		CxVec Result(NewVec,2);
		return Result;
	}
	else
	{
		CxVec Odds = TakeOdds();
		CxVec FTOdds = Odds.fft();
		CxVec Evens = TakeEvens();
		CxVec FTEvens = Evens.fft();
		Complex *FT = new Complex[Size];

		for (unsigned long i = 0; i < Size/2; i++)
		{
			FT[i] = FTEvens.GetEntry(i) + exp(-(2*PI*i/static_cast<double>(Size))*I)*FTOdds.GetEntry(i); 
		}
		for (unsigned long i = Size/2; i < Size; i++)
		{
			FT[i] = FTEvens.GetEntry(i-Size/2) - exp(-(2*PI*(i-Size/2)/static_cast<double>(Size))*I)*FTOdds.GetEntry(i-Size/2); 
		}

		CxVec Result(FT,Size); 
		return Result;
	}
}

