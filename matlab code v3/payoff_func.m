function payoff=payoff_func(X,y)

payoff = max(exp(y)-X, 0);