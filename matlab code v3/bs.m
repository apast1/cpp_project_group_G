Price = 10;
Rate = 0.1;
Volatility = 0.4;
Strike = 15;
Time = 0.5;

Yield=0

[Call, Put] = blsprice(Price, Strike, Rate, Time, Volatility, Yield)