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

	NewDOCallTest("DoCall_BSM_K04.HJH",BSM,4); 
	NewDOCallTest("DoCall_BSM_K06.HJH",BSM,6); 
	NewDOCallTest("DoCall_BSM_K08.HJH",BSM,8); 
	NewDOCallTest("DoCall_BSM_K10.HJH",BSM,10); 
	NewDOCallTest("DoCall_BSM_K12.HJH",BSM,12); 
	NewDOCallTest("DoCall_BSM_K14.HJH",BSM,14); 
	NewDOCallTest("DoCall_BSM_K16.HJH",BSM,16); 
	NewDOCallTest("DoCall_BSM_K20.HJH",BSM,20); 
	NewDOCallTest("DoCall_BSM_K32.HJH",BSM,32); 
	NewDOCallTest("DoCall_BSM_K64.HJH",BSM,64); 
	NewDOCallTest("DoCall_BSM_K128.HJH",BSM,128); 
	NewDOCallTest("DoCall_BSM_K256.HJH",BSM,256); 
	NewDOCallTest("DoCall_BSM_K512.HJH",BSM,512); 
}
