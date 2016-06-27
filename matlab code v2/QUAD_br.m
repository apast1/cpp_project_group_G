function [val,error] =QUAD_br(K,M)



t=0;
dt = (T-t) / M;

k= 2 * (r-Dc)/sig^2 -1;


x = log(S_0/X);


b=zeros(M);

b1=0;





[val,error] =QUAD(S_0,r,sig,X,T,Dc,K)

for i = (M-1):2
    

ymax = log(S_0/X) - 10*sig*sqrt(dt);
ymax = log(S_0/X) + 10*sig*sqrt(dt); % Upper bound of integration range to approximate inf


dy = dt^0.5/K;


N_x= round((ymax-b)/dy) ;



Inv =0;
% for i=1:N_x - 1
%     Inv = Inv + (2*B_func(x,i*dy,sig,dt,k) * payoff_func(X,i*dy) + 4*B_func(x,(i+1/2)*dy,sig,dt,k) * payoff_func(X,(i+1/2)*dy));
% end
% Remain = B_func(x,0,sig,dt,k) * payoff_func(X,0) + B_func(x,dy/2,sig,dt,k) * payoff_func(X,dy/2) +  B_func(x,N_x * dy,sig,dt,k) * payoff_func(X,ymax);

for i=1:N_x - 1
    Inv = Inv + (2*f_func(x,i*dy,sig,dt,k,X) + 4*f_func(x,(i+1/2)*dy,sig,dt,k,X)) ;
end

Remain = f_func(x,0,sig,dt,k,X) + 4 * f_func(x,dy/2,sig,dt,k,X) +  f_func(x,N_x * dy,sig,dt,k,X);

A = exp( (-0.5*k*x-0.125*dt*(sig*k)^2-r*dt) ) / sqrt(2*pi*dt*sig^2);
% A = exp(  -0.5*k*x-0.125*(sig*k)^2*dt - r*dt)  /(2*sig^2*pi*dt)^(1/2);

V = A/6 * dy * (Inv + Remain);

% V = A_func(x,sig,dt,k,r)/6 * dy * (Inv + Remain);
    
val=V;


% Int = f_func(x,0,sig,dt,k,X) + 4* f_func(x,dy/2,sig,dt,k,N_x * dy); %% N_x * dy =? ymax
% 
% 
% for i=1:N_x - 1
% 


% [Call, Put] = blsprice(S_0, X, r, T, sig, Dc);
% 
% error = V - Call;

syms y;


f = exp(-(x-y)^2/(2*sig^2*dt) + 1/2*k*y) * X * (exp(y)-1);
f_4=diff(f,y,4);

% y1=0:1:20;

% y1=exp(1):0.1:20;
y1=0:0.1:ymax;


f_4d=subs(f_4,y,y1);
% plot(y1,f_4d,':b')
max_1= max(f_4d);


error = abs( 1/180 * N_x * dy * (dy)^4 * (  max_1  )   ) *A;

% error = abs( 1/180 * (ymax-b) * (dy/2)^4 * (  max_1  )   )   /8;
