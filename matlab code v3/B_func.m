function B=B_func(x,y,sig,delta_t,k)

B = exp(-(x-y)^2/(2*sig^2*delta_t) + 1/2*k*y);
