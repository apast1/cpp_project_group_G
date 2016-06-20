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
	QuadraticProcess Quadra(10,0.05,0.02,1);

	NewBermudanTest("Bermudan_Qudra_K04.HJH",Quadra,4); 
	NewBermudanTest("Bermudan_Qudra_K06.HJH",Quadra,6); 
	NewBermudanTest("Bermudan_Qudra_K08.HJH",Quadra,8); 
	NewBermudanTest("Bermudan_Qudra_K10.HJH",Quadra,10); 
	NewBermudanTest("Bermudan_Qudra_K12.HJH",Quadra,12); 
	NewBermudanTest("Bermudan_Qudra_K14.HJH",Quadra,14); 
	NewBermudanTest("Bermudan_Qudra_K16.HJH",Quadra,16); 
	NewBermudanTest("Bermudan_Qudra_K20.HJH",Quadra,20); 
	NewBermudanTest("Bermudan_Qudra_K32.HJH",Quadra,32); 
	NewBermudanTest("Bermudan_Qudra_K64.HJH",Quadra,64); 
	NewBermudanTest("Bermudan_Qudra_K128.HJH",Quadra,128); 
	NewBermudanTest("Bermudan_Quadra_K256.HJH",Quadra,256); 
	NewBermudanTest("Bermudan_Quadra_K512.HJH",Quadra,512); 
}
