function [Call] = DOC_MC(S,K,r,sigma,div,T,B,Nchk,NScn)
% Monte Carlo Valuation of Down and Out Call

    dt=T/Nchk; % Time step
    nudt=(r-div-0.5*sigma^2)*dt;
    sigsdt=sigma*sqrt(dt);
    beta = 0.5826;
    B = B*exp(beta*sigma*sqrt(dt)); % BGK(1997,MF) continuity correction
    St=zeros(1,Nchk+1);
    St(1,1) = S;
    Call = 0;
    % Scenario generation
    for i = 1:NScn
        FlagB = 0;
        j=1;
        while (j<= Nchk & ~FlagB)
            St(1,j+1)=St(1,j)*exp(nudt+sigsdt*normrnd(0,1)); % 
            FlagB = (St(1,j+1) < B);
            j = j+1;
        end
        if ~FlagB 
            Call = Call + max(St(1,Nchk+1)-K,0)*exp(-r*T);
        end
    end
    Call = Call / NScn;
