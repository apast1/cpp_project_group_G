clear
clc
S_0 = 100;
r = 0.06;
sig = 0.2;
X = 100;
T = 0.5;
Dc=0;
    
N = 20;      % number of time steps
%     double st_time = clock();

dt = T/N;
% drift = (r - 0.5*sig*sig)*dt;

drift = r * dt;
sgrt = sig*dt^0.5;
discount = exp(-r*T);

acc_vals = 0.0;
acc_squs = 0.0;


K=100000;
D=N;
se_number = 6;


uni_normal=normrnd(0,1,K,N,se_number);


for qq=1:1
    if qq==1
    quasi_normal=uni_normal;
    end

    dd=10;  %M interval number
    for ww=1:dd
        M = K/dd*ww;  % number of sample paths
        acc_se_vals = 0;
        acc_se_squs  = 0;
        for se_i=1:se_number  %samle se
            acc_vals=0;            
            for j=1:M
                S = S_0;
                for i=1:N
                    w=quasi_normal(j,i,se_i);   %%%%
%                   w = normrnd(0,1);

                    S = S*(1 + drift + sgrt*w);
              %     S = S*exp(drift + sgrt*w);
                end
                payoff =discount * max(0, S - X);
                acc_vals = acc_vals + payoff;
            end
%             ccc(se_i,ww) = acc_vals/M;
            c = acc_vals/M;
            acc_se_vals = acc_se_vals + c;
            acc_se_squs = acc_se_squs + c*c;
        end
%         cc(qq,ww)=acc_se_vals/se_number;
%         se(qq,ww)=(acc_se_squs/se_number-(acc_se_vals/se_number)^2)^0.5;
%         se(qq,ww)=std(ccc(:,ww));

    end

end

c
