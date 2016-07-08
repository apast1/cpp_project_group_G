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

I_1=\int_{y_{1}^{min}}^{y_{1}^{max}}\int_{y_{2}^{min}}^{y_{2}=b_m(y_{1})}B(x_1,x_2,y_1,y_2)*\Lambda (y_1,y_2)dy_1dy_2
\newline
I_2=\int_{y_{1}^{min}}^{y_{1}^{max}}\int_{y_{2}=b_m(y_{1})}^{y_{2}^{max}}B(x_1,x_2,y_1,y_2)*\Lambda (y_1,y_2)dy_1dy_2

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


V(x,t)=e^{-r\Delta t} \int_{-\infty }^{\infty}f_{\Delta t}(y|x)V(y,t+\Delta t)dy
\newline

f_{\Delta t}(y|x)=\frac{1}{y \sigma \sqrt{2\pi\Delta t}}exp[-\frac{1}{2}(\frac{ln(y)-ln(x)+(r-\sigma^2/2)\Delta t}{21})]
\newline

p_{x}^{L}(y|x)=\frac{1}{y \sigma \sqrt{2\pi\Delta t}}exp[-\frac{1}{2}(\frac{ln(y)-ln(x)+(r-\sigma^2/2)\Delta t}{21})]




f_{\Delta t}(y|x)=\frac{1}{y \sigma \sqrt{2\pi\Delta t}}exp[-\frac{1}{2}(\frac{ln(y)-ln(x)+(r-\sigma^2/2)\Delta t}{21})]
\newline

f_{\Delta t}(y|x)=P_{X}^{L}(t,y|x)=\frac{P_{Y}^{L}(t,\frac{y^{1-\beta}}{\sigma(1-\beta)}|\frac{x^{1-\beta}}{\sigma(1-\beta)})}{\sigma x^\beta}
\newline

P_{Y}^{L}(t,y|x)=\frac{1}{\sqrt{t}}\phi  ( \frac{y-x}{\sqrt{t}})exp[r(1-\beta)(y^2-x^2)+\frac{\beta}{2(\beta -1)}ln(\frac{y}{x})]\sum_{i=0}^{2}c_{i}(y|x)\frac{t^i}{i!}
\newline

c_{0}(y|x)=1
\newline
c_{1}(y|x)=-\frac{r^2 (1-\beta)^2 (y^3+y^2x+yx^2+x^3)}{6}-\frac{(20r\beta^2-16r\beta-8r\beta^3+\beta^2+4r)(y+x))}{8(\beta-1)^2}-\frac{\beta}{4(\beta-1)^2yx}
f_{\Delta t}(y|x)=P_{X}^{L}(t,y|x)=\frac{P_{Y}^{L}(t,ln\frac{y^{1-\beta}}{\sigma(1-\beta)}|\frac{x^{1-\beta}}{\sigma(1-\beta)})}{\sigma x^\beta}

\newline

f_{\Delta t}(y|x)=P_{X}^{L}(t,y|x)=\frac{P_{Y}^{L}(t,\frac{1}{b\sigma}ln\frac{y}{y+b}|\frac{1}{b\sigma}ln\frac{x}{x+b})}{\sigma x^\beta}
\newline
P_{Y}^{L}(t,y|x)=\frac{1}{\sqrt{t}}\phi  ( \frac{y-x}{\sqrt{t}})\sum_{i=0}^{2}c_{i}(y|x)\frac{t^i}{i!}exp[r(1-\beta)(y^2-x^2)+\frac{\beta}{2(\beta -1)}ln(\frac{y}{x})]

\newline
P_{Y}^{L}(t,y|x)=\frac{1}{\sqrt{t}}\phi  ( \frac{y-x}{\sqrt{t}})\sum_{i=0}^{2}c_{i}(y|x)\frac{t^i}{i!}exp(\frac{(-e^{by\sigma}+e^{bx\sigma})r}{b^2\sigma^2}+\frac{(y-x)(2r-b^2\sigma^2)}{b\sigma}+ln(\frac{1-e^{by\sigma}}{1-e^{bx\sigma}}))
\newline

c_{1}(y|x)=-\frac{r^2 (1-\beta)^2 (y^3+y^2x+yx^2+x^3)}{6}-\frac{(20r\beta^2-16r\beta-8r\beta^3+\beta^2+4r)(y+x))}{8(\beta-1)^2}-\frac{\beta}{4(\beta-1)^2yx}



\newline
P_{Y}^{L}(t,y|x)=\frac{1}{\sqrt{t}}\phi  ( \frac{y-x}{\sqrt{t}})\sum_{i=0}^{2}c_{i}(y|x)\frac{t^i}{i!}exp(\frac{(-e^{by\sigma}+e^{bx\sigma})r}{b^2\sigma^2}+\frac{(y-x)(2r-b^2\sigma^2)}{b\sigma}+ln(\frac{1-e^{by\sigma}}{1-e^{bx\sigma}}))
\newline

c_{1}(y|x)=-\frac{r^2 (1-\beta)^2 (y^3+y^2x+yx^2+x^3)}{6}-\frac{(20r\beta^2-16r\beta-8r\beta^3+\beta^2+4r)(y+x))}{8(\beta-1)^2}-\frac{\beta}{4(\beta-1)^2yx}
\newline

c_{1}(y|x)=-\frac{r^2 (e^{2by\sigma}-e^{2bx\sigma})}{4b^3\sigma^3(y-x)}-\frac{r(r+b^2\sigma^2)(e^{by\sigma}-e^{bx\sigma})}{b^3\sigma^3(y-x)}-\frac{(-2r+b^2\sigma^2)^2}{8b^2\sigma^2}



V_i=X\cdot (e^{y_i}-1)    \quad i=0,1...,2N
\newline
y_i = bm+i \delta y
\newline

\sum_{i=2I_{m,x}^{-}}^{2I_{m,x}^{+}}f_{simpson_+}(V(t_{m+1},x_i))  =\sum_{i=0}^{2I_{m,x}^{+}}f_{simpson}(V(t_{m+1},x_i))+\sum_{i=-2I_{m,x}^{-}}^{0}(X-V(t_{m+1},x_i))


dS_t=rS_tdt+\sqrt{v_t}dW_t^{(1)}
\newline
dv_t=\kappa (\theta -v_t)dt+\sigma _V\sqrt{v_t}dW_t^{(2)}
\newline
dW_t^{(1)}dW_t^{(1)} = \rho dt

\newline
v_t^{max}=v_t + \xi\sqrt{v_tt}
\newline

v_t^{max}=v_t - \xi\sqrt{v_tt}

\newline


x_t^{max}=x_t + (r- v_t/2)t + \xi \sqrt{v_tt}
\newline
x_t^{min}=x_t + (r- v_t/2)t - \xi \sqrt{v_tt}
%	return LogAsset + (rate - 0.5*Var)*Time + Deviations*sqrt(Var)*sqrt(Time);
% 	double Guess = vValue + Deviations*sqrt(vValue*Time);

\newline
