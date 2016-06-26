function val = QUAD_DO(S,E,BAR,sig,r,Dc,T,M)
dt = T/M;
k= 2 * (r-Dc)/sig^2 -1;
BAR(M)=max(BAR(M), E);
dy = sqrt(dt)/4;
for j = M+1:1
    if j>1
        ymax = log(S/E) + 10 sig* sqrt(j*dt);
    end
    
    b(j) = log(bar(j)/E);
    Nplus(j) = around((ymax-b(j))/dy);
end

Nplus(1)=0;

V = zeros(max(Nplus))
for j=M+1:1
    V1 = V
    for i = 1:Nplus(j)
        for ii=[1,2]
            x = b(j)+(i+0.5* (ii-1) )*dy
            if j==M+1
                x = log(S/E)
            else
            
                A = A = exp( (-0.5*k*x-0.125*dt*(sig*(k+2))^2-r*dt) ) / sqrt(2*pi*dt*sig^2);
                qstar = 10*sig*sqrt(dt);
                iplus = min(aroud(   (x + qstar -b(j+1)  ) / dy ), Nplus(j+1));
                iminus = max (aroud ((x - qstar -b(j+1)  ) / dy) , 0);
                int = 0;

                for iii = iplus:-1:iminus

                    for iv = [2,1]

                        y = b(j+1) + (iii + 0.5* (iv-1) ) *dy
                        Bxy = exp(-(x-y)^2/(2*sig^2*dt) + 1/2*k*y)
                        f = Bxy * V1(iii,iv)
                        int = int + (dy/6) * (2+ 2 * iv)*f
                    end
                end

                V(i,ii) = a*(int - (dx/6)*f)
            
            end
            
            if j=1
                break
            end
        end
    end
end

val = V(0,0)
            
            
            
            
            
            
            
            
            
            