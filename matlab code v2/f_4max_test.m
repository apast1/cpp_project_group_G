S_0 = 100;
r = 0.06;
sig = 0.2;
X = 90;
T = 0.5;
Dc=0;




t=0;
dt = T-t;
k= 2 * (r-Dc)/sig^2 -1;
    
x = log(S_0/X);
% ymin = -10;

b = log(X/X);


% ymax = 100*log(S_0);

% dy = sqrt(dt) / 4;


ymax = log(S_0/X) + 10*sig*sqrt(dt); % Upper bound of integration range to approximate inf
% ymax =100;


syms y;


f = exp(-(x-y)^2/(2*sig^2*dt) + 1/2*k*y) * X * (exp(y)-1);
f_4=diff(f,y,4);

% y1=0:1:20;

% y1=exp(1):0.1:20;
y1=0:0.1:ymax;


f_4d=subs(f_4,y,y1);
plot(y1,f_4d,':b')
max_1= max(f_4d);