#include <iostream>

template <typename T>
inline const T& maximum(const T& x,const T& y)
{
   if(y > x){
      return y;
   }
   else{
      return x;
   }
}

int main(void)
{
   using namespace std;
   int a=3,b=7;
   float x=3.0,y=7.0;
   //Calling template function
   std::cout << maximum<int>(a,b) << std::endl;         //?? 7
   std::cout << maximum(x, b) << std::endl;             //????????
   std::cout << maximum<double>(x,y) << std::endl;  //?? 7
   return 0;
}
