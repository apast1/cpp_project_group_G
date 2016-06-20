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
/*
	string name = "sbc_o9";
	name += "_s100b95_dr1";
	name += ".csv";	

	for (unsigned long V = 1; V <= 16; V++)
		{

			NewDOCallTest(name, SABR, V*4*3); 
		}
*/
//bermudan put test

	string name = "sbp";
	name += "_o25_s105m10";
	name += ".csv";	


for (unsigned long V = 1; V <= 16; V++)
{

//	string name1 = "sbp";
//	name1 += "_o10_p";
//	name1 += ".csv";	

	NewBermudanTest(name, SABR, V*4*5,25); 
//	NewBermudanTest(name1, SABR, V*4,10); 

}

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
