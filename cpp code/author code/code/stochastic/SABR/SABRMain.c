#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>
#include <time.h> 
#include <stdio.h>
#include <string.h>


#include "SABR2.h"
#include "TestsB2.h"
#include "SVQUAD2.h"


int main()
{
	SABRProcess SABR(log(100.0),0.2,0.05,0.5,0.1,-0.75);
//barrier call test

for (unsigned long V = 1; V <= 32; V++)
	{
		string name = "sbc_v4o5_intel_p_smp_1";
		name += "_";
		name += "100.csv";	

		NewDOCallTest(name, SABR, V*4); 
	}

//bermudan put test
/*
for (unsigned long V = 2; V < 3; V++)
{
	for (int W = 8; W <=9; W++)
	{
		string name = "sbp";
		name += "_ppr_o10";
		name += ".csv";	

		NewBermudanTest(name, SABR, static_cast<unsigned long>(pow(2.0,W)), V*5); 
	}
}
*/
//euro call test
/*for (unsigned long V = 1; V < 51; V++)
	{
		string name = "ppsec";
		name += ".csv";	

		EuroTest(name,SABR, V*4, 0.1,98); 
	}*/
//euro put test
/*
for (unsigned long V = 1; V < 100; V++)
	{
		string name = "sec_v4_5";
		name += ".csv";	

		EuroPutTest(name,SABR, V*4, 0.1,102); 
	}
*/
}
