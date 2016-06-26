% 
% for i=1:N_x - 1
%     Inv = Inv + (2*f_func(x,i*dy,sig,dt,k,X) + 4*f_func(x,(i+1/2)*dy,sig,dt,k,X)) ;
% end
% 
% Remain = f_func(x,0,sig,dt,k,X) + 4 * f_func(x,dy/2,sig,dt,k,X) +  f_func(x,N_x * dy,sig,dt,k,X)



a =-1;
b= 8;
dx=0.1
N_x=round ((b-a)/dx)

Inv=0
for i=1:N_x - 1
    Inv = Inv + (2*f_func_x(a+i*dx) + 4*f_func_x( a+(i+1/2) *dx)) ;
end

Remain = f_func_x(a) + 4 * f_func_x(a+dx/2) +  f_func_x(N_x*dx)


val = dx/6*(Remain+Inv)

Q = quad('exp(x)',a,b)