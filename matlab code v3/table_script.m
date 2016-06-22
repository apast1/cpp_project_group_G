
S_0=100;
Dc=0;

K=10;


data= [95,0.2,0.06,0.5;105,0.4,0.1,1];
% [X,sigma,r,T]
test_matrix = trans4_16(data);
val_m=zeros(16,1);
error_m=zeros(16,1);

for i=1:16
    X = test_matrix(i,1);
    sig = test_matrix(i,2);
    r = test_matrix(i,3);
    T = test_matrix(i,4);
    [val,error]=QUAD(S_0,r,sig,X,T,Dc,K);
    [Call, Put] = blsprice(S_0, X, r, T, sig, Dc);
    val1 = Call;
    val_m(i)=val;
    error_m(i)=error;
    r(i)=val-val1;
end

RMSE = std(error_m)



