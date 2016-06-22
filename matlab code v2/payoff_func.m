function payoff=payoff_func(X,y)

payoff = X*max(exp(y)-1, 0);