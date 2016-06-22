
T = 0.2;
Dc=0;

K=30;


data= [95,0.2,0.06,0.5;105,0.4,0.1,1];
% [S_0,X,sigma,r]
test_matrix = trans4_16(data);
val_m=zeros(16,1);
error_m=zeros(16,1);

for i=1:16
    S_0 = test_matrix(i,1);
    X = test_matrix(i,2);
    r = test_matrix(i,3);
    sig = test_matrix(i,4);
    [val,error]=QUAD(S_0,r,sig,X,T,Dc,K);
    val_m(i)=val;
    error_m(i)=error;
end

RMSE = std(error_m)
    



