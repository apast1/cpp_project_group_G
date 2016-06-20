#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>
#include <time.h> 

#include "BSProcess.h"
#include "CEV.h"
#include "Tests.h"
#include "BermudanPut.h"
#include "DownAndOutCalls.h"
#include "Quadratic.h"

int main()
{
	BSProcess BSM(10,0.05,0.2);

	NewDOCallTest("docall_const_K04.csv",BSM,4); 
	NewDOCallTest("docall_const_K06.csv",BSM,6); 
	NewDOCallTest("docall_const_K08.csv",BSM,8); 
	NewDOCallTest("docall_const_K10.csv",BSM,10); 
	NewDOCallTest("docall_const_K16.csv",BSM,16); 
	NewDOCallTest("docall_const_K32.csv",BSM,32); 
	NewDOCallTest("docall_const_K64.csv",BSM,64); 
	NewDOCallTest("docall_const_K128.csv",BSM,128); 
}
