#ifndef TESTS_H
#define TESTS_H

#include <fstream>
#include <iomanip>
#include <ctime>

#include "LVProcess.h"
#include "BermudanPut.h"
#include "DownAndOutCalls.h"

using namespace std;

const double DEVS = 7.50;

void BermudanTest(char* Name, LVProcess& Proc, unsigned long N);

void NewBermudanTest(char* Name, LVProcess& Proc, unsigned long K);

void DOCallTest(char* Name, LVProcess& Proc, unsigned long N);

//////////////////////////////////////////////////////////////////////////////

void NewBermudanTest(char* Name, LVProcess& Proc, unsigned long K)
{
	ofstream filu;
	filu.open(Name);
	filu.precision(15); 
	cout.precision(15);
	
	double Maturity;
	unsigned long Obs;
	double Strike;

	clock_t Alku, Loppu;
	vector<double> Times; 
	double Step = 0.0;
	double Price = 0.0;
	double Time = 0.0;
	unsigned long Kierrokset = 100;
	
	if (K > 10)
	{
		Kierrokset = 10;
	}
	if (K > 64)
	{
		Kierrokset = 1;
	}

	Maturity = 0.5;
	Strike = 9.5;
	Obs = 6; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewBermudanPut(Strike,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 0.5;
	Strike = 9.5;
	Obs = 30; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewBermudanPut(Strike,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  

	Maturity = 0.5;
	Strike = 10.5;
	Obs = 6; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewBermudanPut(Strike,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 0.5;
	Strike = 10.5;
	Obs = 30; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewBermudanPut(Strike,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 1.0;
	Strike = 9.5;
	Obs = 6; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewBermudanPut(Strike,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 1.0;
	Strike = 9.5;
	Obs = 30; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewBermudanPut(Strike,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 1.0;
	Strike = 10.5;
	Obs = 6; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewBermudanPut(Strike,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 1.0;
	Strike = 10.5;
	Obs = 30; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewBermudanPut(Strike,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	filu.close(); 
}

void BermudanTest(char* Name, LVProcess& Proc, unsigned long N)
{
	ofstream filu;
	filu.open(Name);
	filu.precision(15); 
	cout.precision(15);
	
	double Maturity;
	unsigned long Obs;
	double Strike;

	clock_t Alku, Loppu;
	vector<double> Times; 
	double TimeStep = 0.0;
	double Price = 0.0;
	double Time = 0.0;

	Maturity = 0.5;
	Strike = 9.5;
	Obs = 6; 
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = BermudanPut(Strike,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 0.5;
	Strike = 9.5;
	Obs = 30; 
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = BermudanPut(Strike,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  

	Maturity = 0.5;
	Strike = 10.5;
	Obs = 6; 
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = BermudanPut(Strike,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 0.5;
	Strike = 10.5;
	Obs = 30; 
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = BermudanPut(Strike,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 1.0;
	Strike = 9.5;
	Obs = 6; 
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = BermudanPut(Strike,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 1.0;
	Strike = 9.5;
	Obs = 30; 
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = BermudanPut(Strike,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 1.0;
	Strike = 10.5;
	Obs = 6; 
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = BermudanPut(Strike,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Maturity = 1.0;
	Strike = 10.5;
	Obs = 30; 
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = BermudanPut(Strike,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << N << "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 
	filu.close(); 
}

void DOCallTest(char* Name, LVProcess& Proc, unsigned long N)
{
	ofstream filu;
	filu.open(Name);
	filu.precision(15); 
	cout.precision(15);

	unsigned long Obs;
	double Strike;
	double Barrier;
	double Maturity = 0.5;

	clock_t Alku, Loppu;
	vector<double> Times; 
	double TimeStep = 0.0;
	double Price = 0.0;
	double Time = 0.0;

	Obs = 6; 
	Strike = 10.0;
	Barrier = 9.0;
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = DOCall(Strike,Barrier,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  
	cout << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  

	Obs = 6; 
	Strike = 10.0;
	Barrier = 9.9;
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = DOCall(Strike,Barrier,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  
	cout << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  

	Obs = 6; 
	Strike = 10.5;
	Barrier = 9.0;
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = DOCall(Strike,Barrier,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  
	cout << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  

	Obs = 6; 
	Strike = 10.5;
	Barrier = 9.9;
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = DOCall(Strike,Barrier,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  
	cout << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  

	Obs = 125; 
	Strike = 10.0;
	Barrier = 9.0;
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = DOCall(Strike,Barrier,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  
	cout << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  

	Obs = 125; 
	Strike = 10.0;
	Barrier = 9.9;
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = DOCall(Strike,Barrier,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  
	cout << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  

	Obs = 125; 
	Strike = 10.5;
	Barrier = 9.0;
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = DOCall(Strike,Barrier,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  
	cout << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  

	Obs = 125; 
	Strike = 10.5;
	Barrier = 9.9;
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	Price = DOCall(Strike,Barrier,Times,Proc,N); 
	Loppu = clock();
	Time = difftime(Loppu,Alku); 
	filu << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  
	cout << N << "," << Obs << "," << Strike << "," << Barrier << "," << Price << "," << Time << "\n";  

	filu.close();
}

void NewDOCallTest(char* Name, LVProcess& Proc, unsigned long K)
{
	ofstream filu;
	filu.open(Name);
	filu.precision(15); 
	cout.precision(15);
	
	double Maturity = 0.5;
	unsigned long Obs;
	double Strike;
	double Barrier;

	clock_t Alku, Loppu;
	vector<double> Times; 
	double Step = 0.0;
	double Price = 0.0;
	double Time = 0.0;
	unsigned long Kierrokset = 100;
	
	if (K > 10)
	{
		Kierrokset = 10;
	}
	if (K > 64)
	{
		Kierrokset = 1;
	}

	Strike = 10.0;
	Barrier = 9.5;
	Obs = 6; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewDOCall(Strike,Barrier,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
 
	Strike = 10.0;
	Barrier = 9.5;
	Obs = 30; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewDOCall(Strike,Barrier,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  

	Strike = 10.0;
	Barrier = 9.9;
	Obs = 6; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewDOCall(Strike,Barrier,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  

	Strike = 10.0;
	Barrier = 9.9;
	Obs = 30; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewDOCall(Strike,Barrier,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  

	Strike = 10.5;
	Barrier = 9.5;
	Obs = 6; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewDOCall(Strike,Barrier,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  

	Strike = 10.5;
	Barrier = 9.5;
	Obs = 30; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewDOCall(Strike,Barrier,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  

	Strike = 10.5;
	Barrier = 9.9;
	Obs = 6; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewDOCall(Strike,Barrier,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  

	Strike = 10.5;
	Barrier = 9.9;
	Obs = 30; 
	Step = 7.5*Maturity/log(static_cast<double>(Obs))/static_cast<double>(K);
	Times = SplitTime(Maturity,Obs); 
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = NewDOCall(Strike,Barrier,Times,Proc,Step); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	filu << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  

	filu.close(); 
}

#endif
