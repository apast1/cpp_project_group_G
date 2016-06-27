Sigma           = 0.1;
AssetPrice      = 100;
OptSpec         = 'call';
Strike          = 105;
Settle          = '01-Jan-2003';
ExerciseDates   = '01-jan-2007';
BarrierSpec     = 'DO'; % up and in
Barrier         = 98;
AmericanOpt     = 0;
DividendType    = []; 
DividendAmounts = 0;
ExDividendDates = []; 

StockSpec       = stockspec(Sigma, AssetPrice, DividendType, ...
                  DividendAmounts, ExDividendDates);
RateSpec        = intenvset('Rates', 0.05, 'StartDates', '01-Jan-2003', ...
                  'EndDates', '01-Jan-2007', 'Compounding', -1);

ValuationDate   = '1-Jan-2003';
Maturity        = '01-Jan-2007';


TimeSpecCRR     = crrtimespec(ValuationDate, Maturity, 100);
CRRTree         = crrtree(StockSpec, RateSpec, TimeSpecCRR);
PriceCRR        = barrierbycrr(CRRTree, OptSpec, Strike, Settle, ...
                  ExerciseDates, AmericanOpt, BarrierSpec, Barrier)


TimeSpecEQP     = eqptimespec(ValuationDate, Maturity, 100);
EQPTree         = eqptree(StockSpec, RateSpec, TimeSpecEQP);
PriceEQP        = barrierbyeqp(EQPTree, OptSpec, Strike, Settle, ...
                  ExerciseDates, AmericanOpt, BarrierSpec, Barrier)
              
              
              
bar = Barrier;
r= 0.05;
T=4;
Dc =0;
M=10;

V1=DOC_MC(AssetPrice,Strike,r,Sigma,Dc,T,bar,10,1000)
val = QUAD_DO(AssetPrice,Strike,bar,Sigma,r,Dc,T,100)
