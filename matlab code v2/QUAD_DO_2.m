% function val = QUAD_DO_2(S,E,bar,sig,r,Dc,T,M)
% val = QUAD_DO(AssetPrice,Strike,bar,Sigma,r,Dc,4,10)
% clear
S0 = 100;
E = 105;
bar = 98;
sig = 0.1;
Dc = 0;
r = 0.05;
T = 4;
M =10;


dt = T/M;
k= 2 * (r-Dc)/sig^2 -1;
bar = ones(1,M+1)*bar;
bar(M+1)=max(bar(M+1), E);
dy = sqrt(dt)/4;
% b=[];
for j = (M+1):-1:1
    if j>1
        ymax = log(S0/E) + 10 * sig* sqrt( (j-1) *dt);
    end
    
    b(j) = log(bar(j)/E);
    Nplus(j) = round((ymax-b(j))/dy);
end

Nplus(1)=0;

V = zeros(2*max(Nplus),1);
% for j=M+1
for j=M+1:-1:1

    V1 = V;
    for i = 1:2*Nplus(j)+1

            z = b(j)+ (i-1) *0.5 *dy;
            
            if j==M+1
                V(i) = E * (exp(z) -1  ) ;
            else
                
                if j==1
                    z = log(S0/E);
                end
                
                A = exp( (-0.5*k*z-0.125*dt*(sig*(k+2))^2-r*dt) ) / sqrt(2*pi*dt*sig^2);
                qstar = 10*sig*sqrt(dt);
                iplus = min(round(   (z + qstar -b(j+1)  ) / dy ), Nplus(j+1));
                iminus = max (round ((z - qstar -b(j+1)  ) / dy) , 1);

                int =0;
                
                for iii = 2* iminus : 2 * iplus;

                        y = b(j+1) + iii * 0.5 *dy;
                        Bxy = exp(-(z-y)^2/(2*sig^2*dt) + 1/2*k*y);
                        f = Bxy * V1(iii);
                        
                        if mod(iii,2)==0
                            int = int + 2 * f;
                        else
                            int = int + 4 * f;
                        end
                        
                        if iii == 2* iminus
                            f1=f;
                        elseif  iii == 2* iplus
                            f2=f;
                        end
                        
                end                
                
                V(i) = dy * A / 6 * (int-f1-f2);
            
            end
            
            if j==1
                break
            end

    end
end

val = V(1);
val
  
            
            
            
            