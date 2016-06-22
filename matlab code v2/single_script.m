S_0 = 100;
r = 0.06;
sig = 0.2;
X = 95;
T = 0.1;
Dc=0;

% 
% S_0 = 100;
% r = 0.1;
% sig = 0.4;
% X = 150;
% T = 0.1;
% Dc=0;


K=100;
[val,error] =QUAD(S_0,r,sig,X,T,Dc,K)