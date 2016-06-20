#ifndef TESTS_H
#define TESTS_H

#include <fstream>
#include <iomanip>
#include <ctime>
#include <string>

#include "SVProcess.h"
#include "SVQUAD2.h"
#include "GenericQUAD4.h"

using namespace std;


void NewBermudanTest(string Name, SVProcess& Proc, unsigned long K, unsigned long O);
void NewDOCallTest(string Name, SVProcess& Proc, unsigned long K);
void EuroTest(string Name,SVProcess& Proc, unsigned long K, double Maturity,double Strike);
void EuroPutTest(string Name,SVProcess& Proc, unsigned long K, double Maturity, double Strike);

//////////////////////////////////////////////////////////////////////////////

void NewBermudanTest(string Name, SVProcess& Proc, unsigned long K, unsigned long O)
{
	ofstream filu;
	filu.open(Name.c_str(),ios::app);
	filu.precision(15); 
	cout.precision(15);
//	filu << CLOCKS_PER_SEC<<"\n";	
	double Maturity;
	unsigned long Obs;
	double Strike;

	clock_t Alku, Loppu;
	vector<double> Times; 
	unsigned long Steps = 0;
	double Price = 0.0;
	double Time = 0.0;
	unsigned long Kierrokset = 10;
	
	if (K > 10)
	{
		Kierrokset = 1;
	}

	Maturity = 1.0;
	Strike = 105;
	Obs = O; 

	Times = SplitTime(Maturity,Obs);
	
	Steps = K;
	
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = BermudanPut(Strike,Times,Proc,Steps); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 

	filu << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Maturity << "," << Strike << "," << Obs << "," << Price << "," << Time << "\n";  
 

	filu.close(); 
}

void NewDOCallTest(string Name, SVProcess& Proc, unsigned long K)
{
	ofstream filu;
	filu.open(Name.c_str(),ios::app);
	filu.precision(15); 
	cout.precision(15);
//	filu << CLOCKS_PER_SEC<<"\n"; 	
	double Maturity = 0.5;
	unsigned long Obs;
	double Strike;
	double Barrier;

	clock_t Alku, Loppu;
	vector<double> Times; 
	unsigned long Steps = 0;
	double Price = 0.0;
	double Time = 0.0;
	unsigned long Kierrokset = 10;
	
	if (K > 10)
	{
		Kierrokset = 1;
	}

	Strike = 100;
	Barrier = 95;
	Obs = 9; 

	Times = SplitTime(Maturity,Obs);
	Steps = K;
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = DOCall(Strike,Barrier,Times,Proc,Steps); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset);


	filu << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  
	cout << K <<  "," << Strike << "," << Barrier << "," << Obs << "," << Price << "," << Time << "\n";  

	filu.close(); 
}

void EuroTest(string Name,SVProcess& Proc, unsigned long K, double Maturity,double Strike)
{
	ofstream filu;
	filu.open(Name.c_str(),ios::app);
	filu.precision(15);
	cout.precision(15);
//	filu << CLOCKS_PER_SEC<<"\n";
	clock_t Alku, Loppu;
	double Times = 0.0; 
	unsigned long Steps = 0;
	double Price = 0.0;
	double Time = 0.0;
	unsigned long Kierrokset = 10;
	
	if (K > 10)
	{
		Kierrokset = 1;
	}


	Times = Maturity; 
	Steps = K;
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = EuroCall(Strike,Times,Proc,Steps); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
	

	filu << K <<  "," << Maturity << "," << Strike << "," << Price << "," << Time << "\n";
	cout << K <<  "," << Maturity << "," << Strike << "," << Price << "," << Time << "\n";  
 
	filu.close();
}

void EuroPutTest(string Name, SVProcess& Proc, unsigned long K, double Maturity, double Strike)
{

	ofstream filu;
	filu.open(Name.c_str(),ios::app);
	filu.precision(15);
	cout.precision(15);
//	filu << CLOCKS_PER_SEC<<"\n";	
	clock_t Alku, Loppu;
	double Times = 0.0; 	unsigned long Steps = 0;
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

	Times = Maturity;
	Steps = K;
	Alku = clock();
	for (unsigned long V = 0; V < Kierrokset; V++)
	{
		Price = EuroPut(Strike,Times,Proc,Steps); 
	}
	Loppu = clock();
	Time = difftime(Loppu,Alku)/static_cast<double>(Kierrokset); 
 

	filu << K <<  "," << Maturity << "," << Strike << "," << Price << "," << Time << "\n";
	cout << K <<  "," << Maturity << "," << Strike << "," << Price << "," << Time << "\n";  

	filu.close();

}
#endif
