% function val = QUAD_DO_2(S0,E,bar,sig,r,Dc,T,M,K)
% val = QUAD_DO(AssetPrice,Strike,bar,Sigma,r,Dc,4,10)
clear


S0 = 100;
E = 95;
bar = 70;
sig = 0.2;
Dc = 0;
r = 0.06;
T = 0.5;
K = 4;
M=50;

dt = T/M;
k= 2 * (r-Dc)/sig^2 -1;
bar = ones(1,M+1)*bar;
bar(M+1)=max(bar(M+1), E);
dy = sqrt(dt)/K;
% b=[];
for j = (M+1):-1:1
    if j>1
        ymax = log(S0/E) + 10 * sig* sqrt( (j-1) *dt);
    end
    
    b(j) = log(bar(j)/E);
    Nplus(j) = round((ymax-b(j))/dy);
end

Nplus(1)=0;
Vr=zeros(2*max(Nplus)+1,M+1);
xx=zeros(2*max(Nplus)+1,M+1);

V = zeros(2*max(Nplus)+1,1);
x = zeros(2*max(Nplus)+1,1);
% j=M+1
for j=M+1:-1:1
%     for j=M+1:2

    V1 = V;
    V=[];
    x=[];
    for i = 1:2*Nplus(j)+1

            x(i) = b(j)+ (i-1) *0.5 *dy;
            
            if j==M+1
                V(i) = E * (exp(x(i)) -1  ) ;
            else
                
                if j==1
                    x(i) = log(S0/E);
                end
                
                A = exp( (-0.5*k*x(i)-0.125*dt*(sig*(k+2))^2-r*dt) ) / sqrt(2*pi*dt*sig^2);
                
                qstar = 10*sig*sqrt(dt);
                iplus = min(round(   (x(i) + qstar -b(j+1)  ) / dy ), Nplus(j+1));
%                 iplus = max(min(round(   (x(i) + qstar -b(j+1)  ) / dy ), Nplus(j+1)),0);
                iminus = max (round ((x(i) - qstar -b(j+1)  ) / dy) , 0);
%                 iplus = 5;
%                 iminus=3;
                int =0;
                
                
                
                for iii = 2* iminus : 2 * iplus;
%                                     for iii =  2 * iplus : -1:2* iminus;
%                                     for iii =  2 * iplus;

                        
                        y = b(j+1) + iii * 0.5 *dy;
                        Bxy = exp(-(x(i)-y)^2/(2*sig^2*dt) + 1/2*k*y);
                        
                        f = Bxy * V1(iii+1);
                        
                        if mod(iii,2)==0
                            int = int + 2 * f;
                        else
                            int = int + 4 * f;
                        end
%                         f1=0;
                        if iii == 2* iminus
                            f1=f;
                        end
                        
                        if  iii == 2* iplus
                            f2=f;
                        end
                        
                end                
                
                V(i) = dy * A / 6 * (int-f1-f2);
            
            end
            
%             if j==1
%                 break
%             end

    end
    
%     Vr(1:length(V),j)=V;
%     xx(1:length(x),j)=x;
    
end

val = V(1);
val
  
            
            
            
            