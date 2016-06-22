function A=A_func(x,sig,delta_t,k,r)

A = 1/(2*sig^2*pi*delta_t)^1/2 * exp(-1/2*k*x - 1/8*sig^2*(k+2)^2*delta_t - r*delta_t);

