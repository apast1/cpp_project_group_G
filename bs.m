Price = 100;
Rate = 0.1;
Volatility = 0.4;
Strike = 150;
Time = 0.25;

Yield=0

[Call, Put] = blsprice(Price, Strike, Rate, Time, Volatility, Yield)