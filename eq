V(S,t)=A(x) \int_{-\infty }^{\infty}B(x,y)V(y,t+\Delta t)dy
\newline
A(x) = \frac{1}{\sigma \sqrt{2\pi\Delta t }}exp(-\frac{kx}{2} - \frac{(\sigma k)^{2}\Delta t}{8}-r\Delta t)
\newline
B(x,y)=exp(-\frac{(x-y)^{2}}{2\sigma^{2}\Delta t}+\frac{ky}{2})
\newline
k = \frac{2(r-D_c)}{\sigma^{2}}-1
\newline
f(x,y)=B(x,y)\Lambda (y,X)
\newline
V(S,t)=A(x) \int_{-\infty }^{\infty}f(x,y)dy
\newline

L(f)=\int_{a}^{b}f(y)dy
\newline
\widetilde{L(f)}=\frac{\delta y}{6}(f(a) +2\sum_{j=1}^{N-1}f(y_{2j}) + 4\sum_{j=0}^{N-1}f(y_{2j+1})+f(b))
\newline
\delta y = \frac{b-a}{N}
\newline
y_{i}=a+\frac{\delta y}{2}\cdot i     \qquad  i=0,1,...,N
\newline
E_{K,2}(f)=L(f)-\widetilde{L(f)}=







L(f)=\int_{a}^{b}f(y)dy
\newline
\widetilde{L(f)}=\frac{\delta y}{6}(f(a) +2\sum_{j=1}^{N-1}f(y_{2j}) + 4\sum_{j=0}^{N-1}f(y_{2j+1})+f(b))
\newline
\delta y = \frac{b-a}{N}
\newline
y_{i}=a+\frac{\delta y}{2}\cdot i     \qquad  i=0,1,...,N
\newline
E_{K,2}(f)=L(f)-\widetilde{L(f)}=-\frac{b-a}{180}(\frac{H}{2})^{4}\cdot f^{(4)}(\xi )  \qquad \xi\in (a,b)

\newline



L(f)=\int_{a}^{b}f(y)dy
\newline
\widetilde{L(f)}=\frac{\delta y}{6}(f(a) +2\sum_{j=1}^{N-1}f(y_{2j}) + 4\sum_{j=0}^{N-1}f(y_{2j+1})+f(b))
\newline
\delta y = \frac{b-a}{N}
\newline
y_{i}=a+\frac{\delta y}{2}\cdot i     \qquad  i=0,1,...,2N
\newline
E_{N,2}(f)=L(f)-\widetilde{L(f)}=-\frac{b-a}{180}(\frac{H}{2})^{4}\cdot f^{(4)}(\xi )  \qquad \xi\in (a,b)
\newline
T_{m}=t+m\Delta t \qquad m=0,1,...,M
\newline 
\Delta t = \frac{T-t}{M}


y_{m,max}= x_{0}+10\sigma\sqrt{T_m-t}
\newline 
y_{m,min} = x_{0}-10\sigma\sqrt{T_m-t}
\newline 
\delta y = \sqrt{\Delta t}/K  \qquad where \quad \Delta t = T-t


\newline 
N_{m}^{+} = \frac{y_{m,max}-b_{m}}{dy/2}
\newline 
N_{m}^{-} = \frac{b_{m}-y_{m,min}}{dy/2}
\newline 

x_{m,j} = b_{m} + \frac{\delta y}{2} \cdot j  \qquad \quad j\in[N^-,N^+]

\newline 
\widehat{I_{m,j}^{+}}   = \frac{x_{j}+10\sigma\sqrt{\Delta t_{m}}-b_{m}}{dy/2}
\newline 
\widehat{I_{m,j}^{-}} = \frac{x_{j}-10\sigma\sqrt{\Delta t_{m}}-b_{m}}{dy/2}
\newline 
I_{m,j}^{+} = min(\widehat{I_{m,j}^{+}}, N_{m+1}^+)
\newline 
I_{m,j}^{-} = max(\widehat{I_{m,j}^{-}}, N_{m+1}^-)
\newline 

y_{i}^{m,j} = b_{m+1}+ \frac{\delta y}{2} \cdot i \qquad \quad i\in[I_{m,j}^{-},I_{m,j}^{+}]


\newline 



y_{i}^{m,j} = b_{m+1}+ \frac{\delta y}{2} \cdot i \qquad \quad i\in[I_{m,j}^{-},I_{m,j}^{+}]
\newline 

V(t_{m},x_j)=max(\widetilde{V(t_{m},x_j)},Payoff(t_{m},x_j))
\newline 
\widetilde{V(t_{m},x_j)}=\sum_{i=I_{m,j}^{-}}^{I_{m,j}^{+}}f_{simpson}(V(t_{m+1},x_i))
\newline 
g(x) = \widetilde{V(t_{m},x)}-Payoff(t_{m},x)
\newline 
b_{m,new}=b_{m,old}+\Delta b_{m,new}

\newline 

y_{m,max} y_{m}^{max}
y_{m,min} y_{m}^{min}

\newline 

\Delta b_{m,new} = \frac{-\Delta b_{m,old}g(b_{m,new})}{g(b_{m,new})-g(b_{m,old})}

