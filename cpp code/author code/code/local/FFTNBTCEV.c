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
	CEVProcess CEV(10,0.5,0.05,0.6);

	NewBermudanTest("Bermudan_CEV_K04.HJH",CEV,4); 
	NewBermudanTest("Bermudan_CEV_K06.HJH",CEV,6); 
	NewBermudanTest("Bermudan_CEV_K08.HJH",CEV,8); 
	NewBermudanTest("Bermudan_CEV_K10.HJH",CEV,10); 
	NewBermudanTest("Bermudan_CEV_K12.HJH",CEV,12); 
	NewBermudanTest("Bermudan_CEV_K14.HJH",CEV,14); 
	NewBermudanTest("Bermudan_CEV_K16.HJH",CEV,16); 
	NewBermudanTest("Bermudan_CEV_K20.HJH",CEV,20); 
	NewBermudanTest("Bermudan_CEV_K32.HJH",CEV,32); 
	NewBermudanTest("Bermudan_CEV_K64.HJH",CEV,64); 
	NewBermudanTest("Bermudan_CEV_K128.HJH",CEV,128); 
	NewBermudanTest("Bermudan_CEV_K256.HJH",CEV,256); 
	NewBermudanTest("Bermudan_CEV_K512.HJH",CEV,512); 
}
