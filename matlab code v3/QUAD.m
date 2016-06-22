function [val,error] =QUAD_ext(S_0,r,sig,X,T,Dc,K)


t=0;
delta_t = T-t;
k= 2 * (r-Dc)/sig^2 -1;
    
x = log(S_0);
% ymin = -10;

b = log(X);
ymax = 100*log(S_0);

delta_y = delta_t^0.5/K;
N_x= floor((ymax-b)/delta_y);

Inv =0;
for i=1:N_x - 1
    Inv = Inv + (2*B_func(x,i*delta_y,sig,delta_t,k) * payoff_func(X,i*delta_y) + 4*B_func(x,(i+1/2)*delta_y,sig,delta_t,k) * payoff_func(X,(i+1/2)*delta_y));
end
Remain = B_func(x,0,sig,delta_t,k) * payoff_func(X,0) + 4*B_func(x,delta_y/2,sig,delta_t,k) * payoff_func(X,delta_y/2) +  B_func(x,N_x * delta_y,sig,delta_t,k) * payoff_func(X,N_x * delta_y);
V = A_func(x,sig,delta_t,k,r)/6 * delta_y * (Inv + Remain);

val=V;


% [Call, Put] = blsprice(S_0, X, r, T, sig, Dc);
% 
% error = V - Call;

syms y;


% f = exp(-(x-y)^2/(2*sig^2*delta_t) + 1/2*k*y) * X * max(exp(y)-1, 0);
f = exp(-(x-y)^2/(2*sig^2*delta_t) + 1/2*k*y) * X * (exp(y)-1);
f_4=diff(f,y,4);

% y1=0:1:20;

% y1=exp(1):0.1:20;
y1=0:0.1:ymax;


f_4d=subs(f_4,y,y1);
% plot(y1,f_4d,':b')
max_1= max(f_4d);




error = abs( 1/180 * N_x * delta_y * (delta_y/2)^4 * (  max_1  )   );
