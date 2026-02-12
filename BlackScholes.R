#Black Scholes Formula and Replications
#S0 - Stock Price, X- Exercise Price, r- Risk free rate, T- Time of expiration, sigma- sd of log returns
# Call option for European: C0 = S0*N(d1) - Xe^(-rT)N(d2)
#d1 = [ln(S0/X)+(r+(\sigma²/2)]/\sigma\root T
#d2 = [ln(S0/X)+(r-(\sigma²/2)]/\sigma\root T

BlackScholes_base <- function(S, X, r, T, sig){
  d1 <- (log(S/X) + (r + sig^2/2)*T) / (sig * sqrt(T))
  d2 <- (log(S/X) + (r - sig^2/2)*T) / (sig * sqrt(T))
  cbs <- S * pnorm(d1) - X * exp(-r*T) * pnorm(d2)
  return(cbs)
}

call <- BlackScholes_base(S=102, X=100, r=0.02, T=0.5, sig=0.2)
call
