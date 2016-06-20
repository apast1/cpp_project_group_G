#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>
#include <time.h> 
#include <stdio.h>
#include <string.h>


#include "Heston.h"
#include "TestsB2.h"
#include "SVQUAD2.h"


int main()
{
	HestonProcess Heston(log(100.0),0.2,0.05,2,0.3,0.6,-0.75);

	cout << CLOCKS_PER_SEC << endl;
// test for error convergency in european call
/*
for (unsigned long V = 2; V <= 50; V++)
	{	
		string name = "T2E2variance";
		name += "004ht";
		name += ".csv";	
		
		EuroTest(name,Heston, V*4, 0.1,98);
	}
*/
// test for maturity impacts on accuracy in european call
/*
for (unsigned long M = 1; M <=50; M++)
	{
		EuroTest("100.csv",Heston,100,0.01*M,9.8);
		EuroTest("accu.csv",Heston, 500, 0.01*M,9.8);
	}	
*/
// test for error convergency in european put
/*
for (unsigned long V = 2; V <= 50; V++)
	{	
		string name = "epc";
		name += "004";
		name += ".csv";	
		
		EuroPutTest(name,Heston, V*4, 0.1,105);
	} 
*/
//Bermudan call convergence towards American/Euro call price
/*
for (unsigned long V = 1; V<=30; V++)
	{
		string name = "BermudanC";
		name += "v128";
		name += ".csv";

		NewBermCallTest(name,Heston,128,V*2);
	}
*/
// Barrier call error convergency
/*
for (unsigned long V = 1; V <=16 ; V++)
	{	
		string name = "Hdc";
		name += "o9_s100b95r";
		name += ".csv";	
		
		NewDOCallTest(name,Heston, V*4*3);
	}

*/
// Bermudan put error convergency

for (unsigned long V = 1; V <= 40; V++)
	{	
		string name = "hbp";
		name += "o25s105m10";
		name += ".csv";	
		
		NewBermudanTest(name,Heston,V*8,25);
	}

// Bermudan put convergence towards american
/*
for (unsigned long V = 30; V<=36; V++)
	{
		string name = "Bermudanput";
		name += "m01q128ref";
		name += ".csv";
		NewBermudanTest(name,Heston,128,V*2);
		NewBermudanTest(name,Heston,192,V*2);
	}
*/
}

