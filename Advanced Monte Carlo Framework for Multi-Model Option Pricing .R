#Monte Carlo Option Pricing
#Pricing of European Options on Dividend Paying stock
#ST=S0e^(r−σ^2/2)T+σWT where W_T~N(0,1)
#Payoff of the call option is max(S_T-K,0) and for put option is max(K-S_T)

# call put option monte carlo
call_put_mc<-function(nSim=100000, tau, r, sigma, S0, K) {
  
  Z <- rnorm(nSim, mean=0, sd=1)
  WT <- sqrt(tau) * Z
  ST = S0*exp((r - 0.5*sigma^2)*tau + sigma*WT)
  
  # price and standard error of call option
  simulated_call_payoffs <- exp(-r*tau)*pmax(ST-K,0)
  price_call <- mean(simulated_call_payoffs)
  sterr_call <- sd(simulated_call_payoffs)/sqrt(nSim)
  # price and standard error of put option
  simulated_put_payoffs <- exp(-r*tau)*pmax(K-ST,0)
  price_put <- mean(simulated_put_payoffs)
  sterr_put <- sd(simulated_put_payoffs)/sqrt(nSim)
  
  output<-list(price_call=price_call, sterr_call=sterr_call, 
               price_put=price_put, sterr_put=sterr_put)
  return(output)
  
}
set.seed(1)
results<-call_put_mc(n=100000, tau=0.5, r=0.02, sigma=0.2, S0=102, K=100)
results
