Sigma           = 0.2;
AssetPrice      = 100;
OptSpec         = 'call';
Strike          = 95;
Settle          = '01-July-2003';
ExerciseDates   = '01-jan-2004';
BarrierSpec     = 'DO'; % up and in
Barrier         = 70;
AmericanOpt     = 0;
DividendType    = []; 
DividendAmounts = 0;
ExDividendDates = []; 

StockSpec       = stockspec(Sigma, AssetPrice, DividendType, ...
                  DividendAmounts, ExDividendDates);
RateSpec        = intenvset('Rates', 0.06, 'StartDates', '01-July-2003', ...
                  'EndDates', '01-Jan-2004', 'Compounding', -1);

ValuationDate   = '1-July-2003';
Maturity        = '01-Jan-2004';


TimeSpecCRR     = crrtimespec(ValuationDate, Maturity, 100);
CRRTree         = crrtree(StockSpec, RateSpec, TimeSpecCRR);
PriceCRR        = barrierbycrr(CRRTree, OptSpec, Strike, Settle, ...
                  ExerciseDates, AmericanOpt, BarrierSpec, Barrier)


TimeSpecEQP     = eqptimespec(ValuationDate, Maturity, 100);
EQPTree         = eqptree(StockSpec, RateSpec, TimeSpecEQP);
PriceEQP        = barrierbyeqp(EQPTree, OptSpec, Strike, Settle, ...
                  ExerciseDates, AmericanOpt, BarrierSpec, Barrier)
              
              
              
bar = Barrier;
r= 0.06;
T=0.5;
Dc =0;
M=30;
K=100;

V1=DOC_MC(AssetPrice,Strike,r,Sigma,Dc,T,bar,M,1000)
val = QUAD_DO_2_dubug(AssetPrice,Strike,bar,Sigma,r,Dc,T,M,K)
% val1 = QUAD_DO(AssetPrice,Strike,bar,Sigma,r,Dc,T,M,K)
