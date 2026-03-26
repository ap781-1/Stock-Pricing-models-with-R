#Black Scholes Formula and Replications
#S0 - Stock Price, K- Exercise Price, r- Risk free rate, T- Time of expiration, sigma- sd of log returns
# Call option for European: C0 = S0*N(d1) - Xe^(-rT)N(d2)
#d1 = [ln(S0/X)+(r+(\sigma²/2)]/\sigma\root T
#d2 = [ln(S0/X)+(r-(\sigma²/2)]/\sigma\root T

# Black-Scholes formula 
black_scholes <- function(S0, K, T, r, sigma) {
  d1 <- (log(S0/K) + (r + sigma^2/2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
    price <- S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  return(price)
}


call <- black_scholes(S0=102, K=100, r=0.02, T=0.5, sigma=0.2)
call

