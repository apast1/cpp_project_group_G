function val = QUAD_DO(S,E,bar,sig,r,Dc,T,M,K)
% val = QUAD_DO(AssetPrice,Strike,bar,Sigma,r,Dc,4,10)

% S = 100;
% E = 105;
% bar = 98;
% sig = 0.1;
% Dc = 0;
% r = 0.05;
% T = 4;
% M =2;


dt = T/M;
k= 2 * (r-Dc)/sig^2 -1;
bar = ones(1,M+1)*bar;
bar(end)=max(bar(end), E);
dy = sqrt(dt)/K;
b=[];
for j = (M+1):-1:1
    if j>1
        ymax = log(S/E) + 10 * sig* sqrt( (j-1) *dt);
    end
    
    b(j) = log(bar(j)/E);
    Nplus(j) = round((ymax-b(j))/dy);
end

Nplus(1)=0;

V = zeros(max(Nplus),2);
for j=M+1:-1:1
% j=M+1

    V1 = V;
    for i = 1:Nplus(j)+1
        for ii=[1,2]
            x = b(j)+( (i-1) +0.5* (ii-1) )*dy;
            if j==M+1
                V(i,ii) = E * (exp(x) -1  ) ;
                
            else
                if j==1
                    x = log(S/E);
                end
                A = exp( (-0.5*k*x-0.125*dt*(sig*(k+2))^2-r*dt) ) / sqrt(2*pi*dt*sig^2);
                qstar = 10*sig*sqrt(dt);
                iplus = min(round(   (x + qstar -b(j+1)  ) / dy ), Nplus(j+1));
                iminus = max (round ((x - qstar -b(j+1)  ) / dy) , 1);
                int = 0;

                for iii = iplus:-1:iminus

                    for iv = [2,1]

                        y = b(j+1) + (iii + 0.5* (iv-1) ) *dy;
                        Bxy = exp(-(x-y)^2/(2*sig^2*dt) + 1/2*k*y);
                        f = Bxy * V1(iii,iv);
                        int = int + (dy/6) * (2+ 2 * (iv-1))*f;
                    end
                end

                V(i,ii) = A*(int - (dy/6)*f);
            
            end
            
            if j==1
                break
            end
        end
    end
end

val = V(1,1);
            
            
            
val
            
            
            
            
            
            