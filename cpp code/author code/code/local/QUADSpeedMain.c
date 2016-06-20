#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>
#include <time.h> 

//#include "FFT.h"
#include "Processes.h"
//#include "MiscMaths.h"
//#include "Multilevel.h"
#include "Option.h"
#include "QUAD.h"
//#include "NewFFT.h"

using namespace std;

//typedef complex<double> Complex; 

int main()
{
	double AccValue = 0;
	double Value = 0;
	double TotalTime = 0;
	double Error = 0; 

for (unsigned long k = 3; k < 11; k++)
{
	unsigned long K = 2*k; 

	unsigned long N = ceil(15*0.2*sqrt(6.0)*K); 

	BSProcess BS1(100,0.05,0.2); 

	vector<double> Times = SplitTime(0.5,6); 
	
	// 6-95-0.2-0.5
	AccValue = BermudanPut(95,Times,BS1,1000);
	time_t Alku = clock(); 
	Value = BermudanPut(95,Times,BS1,N);
	time_t Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	// 6-105-0.2-0.5
	AccValue = BermudanPut(105,Times,BS1,1000);
	Alku = clock(); 
	Value = BermudanPut(105,Times,BS1,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	Times = SplitTime(1,6); 

	// 6-95-0.2-1
	AccValue = BermudanPut(95,Times,BS1,1000);
	Alku = clock(); 
	Value = BermudanPut(95,Times,BS1,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	//  6-105-0.2-1
	AccValue = BermudanPut(105,Times,BS1,1000);
	Alku = clock(); 
	Value = BermudanPut(105,Times,BS1,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	Times = SplitTime(0.5,125);

	N = ceil(15*0.2*sqrt(125.0)*K); 
	// 125-95-0.2-0.5
	AccValue = BermudanPut(95,Times,BS1,1000);
	Alku = clock(); 
	Value = BermudanPut(95,Times,BS1,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	// 125-105-0.2-0.5
	AccValue = BermudanPut(105,Times,BS1,1000);
	Alku = clock(); 
	Value = BermudanPut(105,Times,BS1,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	Times = SplitTime(1,125);

	// 125-95-0.2-1
	AccValue = BermudanPut(95,Times,BS1,1000);
	Alku = clock(); 
	Value = BermudanPut(95,Times,BS1,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	// 125-105-0.2-1
	AccValue = BermudanPut(105,Times,BS1,1000);
	Alku = clock(); 
	Value = BermudanPut(105,Times,BS1,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue);

	BSProcess BS2(100,0.05,0.4); 

	Times = SplitTime(0.5,6); 
	N = ceil(15*0.4*sqrt(6.0)*K); 
	// 6-95-0.4-0.5
	AccValue = BermudanPut(95,Times,BS2,1000);
	Alku = clock(); 
	Value = BermudanPut(95,Times,BS2,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	// 6-105-0.4-0.5
	AccValue = BermudanPut(105,Times,BS2,1000);
	Alku = clock(); 
	Value = BermudanPut(105,Times,BS2,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	Times = SplitTime(1,6); 

	// 6-95-0.4-1
	AccValue = BermudanPut(95,Times,BS2,1000);
	Alku = clock(); 
	Value = BermudanPut(95,Times,BS2,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	//  6-105-0.4-1
	AccValue = BermudanPut(105,Times,BS2,1000);
	Alku = clock(); 
	Value = BermudanPut(105,Times,BS2,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	Times = SplitTime(0.5,125);

	N = ceil(15*0.4*sqrt(125.0)*K); 

	// 125-95-0.4-0.5
	AccValue = BermudanPut(95,Times,BS2,1000);
	Alku = clock(); 
	Value = BermudanPut(95,Times,BS2,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	// 125-105-0.4-0.5
	AccValue = BermudanPut(105,Times,BS2,1000);
	Alku = clock(); 
	Value = BermudanPut(105,Times,BS2,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	Times = SplitTime(1,125); 

	// 125-95-0.4-1
	AccValue = BermudanPut(95,Times,BS2,1000);
	Alku = clock(); 
	Value = BermudanPut(95,Times,BS2,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue); 

	// 125-105-0.4-1
	AccValue = BermudanPut(105,Times,BS2,1000);
	Alku = clock(); 
	Value = BermudanPut(105,Times,BS2,N);
	Loppu = clock(); 
	TotalTime += difftime(Loppu,Alku)/CLOCKS_PER_SEC; 
	Error = (Value-AccValue)*(Value-AccValue);

	cout << "K: " << K << ". RMS: " << sqrt(Error/16.0) << ". Time: " << TotalTime/16.0 << ".\n";  
}
}