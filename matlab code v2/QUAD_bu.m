% function val = QUAD_bu(S0,E,bar,sig,r,Dc,T,M,K)

clear
S0 = 100;
E = 95;
bar = 70;
sig = 0.2;
Dc = 0;
r = 0.06;
T = 0.5;
K = 4;
M=2;

dt = T/M;
k= 2 * (r-Dc)/sig^2 -1;

dy = sqrt(dt)/K;
b=[];
b(M+1)=0;


% Nplus(M+1) = round((ymax-b(M+1))/dy);
for j = (M+1):-1:1
    if j>1
        ymax = log(S0/E) + 10 * sig* sqrt( (j-1) *dt);
        ymin = log(S0/E) - 10 * sig* sqrt( (j-1) *dt);
            flag=j;

    end
    Nplus(j) = round((ymax-ymin)/dy);
end
Nplus(1) = 0;


% %%%%%%%%%%%%% 
% Nplus(1)=0;
% Vr=zeros(2*max(Nplus)+1,M+1);
% xx=zeros(2*max(Nplus)+1,M+1);

V = zeros(2*max(Nplus)+1,1);
x = zeros(2*max(Nplus)+1,1);
% for j=M+1:-1:1
    for j=M+1:1


    V1 = V;
    V=[];
    x=[];
    for i = 1:2*Nplus(j)+1

%             x(i) = b(j)+ (i-1) *0.5 *dy;
            x(i) = log(S0/E) - 10 * sig* sqrt( (j-1) *dt) + (i-1) *0.5 *dy;
            
            
            
            if j==M+1
                x(i) = b(j)+ (i-1) *0.5 *dy;
                V(i) = E * (exp(x(i)) -1  ) ;
            else
                
                if j==1
                    x(i) = log(S0/E);
                end
                
                A = exp( (-0.5*k*x(i)-0.125*dt*(sig*(k+2))^2-r*dt) ) / sqrt(2*pi*dt*sig^2);
                
                qstar = 10*sig*sqrt(dt);
                iplus = min(round(   (x(i) + qstar -b(j+1)  ) / dy ), Nplus(j+1));
                iminus = max (round ((x(i) - qstar -b(j+1)  ) / dy) , 0);

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
  
            
            
            
            