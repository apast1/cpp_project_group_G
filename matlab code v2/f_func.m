function f=f_func(x,y,sig,delta_t,k,X)

f = exp(-(x-y)^2/(2*sig^2*delta_t) + 1/2*k*y) * X*max(exp(y)-1, 0);
