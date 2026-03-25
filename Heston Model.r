# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvtnorm)
library(RColorBrewer)

# Parameters
# simulation dependent
S0 <- 100.0             # asset price
T <- 1.0                # time in years
r <- 0.02               # risk-free rate
N <- 252                # number of time steps in simulation
M <- 1000               # number of simulations

# Heston dependent parameters
kappa <- 3              # rate of mean reversion of variance under risk-neutral dynamics
theta <- 0.20^2         # long-term mean of variance under risk-neutral dynamics
v0 <- 0.25^2            # initial variance under risk-neutral dynamics
rho <- 0.7              # correlation between returns and variances under risk-neutral dynamics
sigma <- 0.6            # volatility of volatility

# Heston model simulation function
heston_model_sim <- function(S0, v0, rho, kappa, theta, sigma, T, N, M) {
  """
  Inputs:
   - S0, v0: initial parameters for asset and variance
   - rho   : correlation between asset returns and variance
   - kappa : rate of mean reversion in variance process
   - theta : long-term mean of variance process
   - sigma : vol of vol / volatility of variance process
   - T     : time of simulation
   - N     : number of time steps
   - M     : number of scenarios / simulations
  
  Outputs:
  - asset prices over time (matrix)
  - variance over time (matrix)
  """
  dt <- T / N
  
  # Arrays for storing prices and variances
  S <- matrix(S0, nrow = N + 1, ncol = M)
  v <- matrix(v0, nrow = N + 1, ncol = M)
  
  # Correlation matrix
  cov_matrix <- matrix(c(1, rho, rho, 1), nrow = 2)
  
  # Sampling correlated brownian motions under risk-neutral measure
  for (i in 1:N) {
    # Generate correlated random numbers
    Z <- rmvnorm(M, mean = c(0, 0), sigma = cov_matrix)
    
    # Update asset prices and variance
    S[i + 1, ] <- S[i, ] * exp((r - 0.5 * v[i, ]) * dt + 
                                  sqrt(v[i, ] * dt) * Z[, 1])
    v[i + 1, ] <- pmax(v[i, ] + kappa * (theta - v[i, ]) * dt + 
                          sigma * sqrt(v[i, ] * dt) * Z[, 2], 0)
  }
  
  return(list(S = S, v = v))
}

# Simulate with different correlations
rho_p <- 0.98
rho_n <- -0.98

result_p <- heston_model_sim(S0, v0, rho_p, kappa, theta, sigma, T, N, M)
result_n <- heston_model_sim(S0, v0, rho_n, kappa, theta, sigma, T, N, M)

S_p <- result_p$S
v_p <- result_p$v
S_n <- result_n$S
v_n <- result_n$v

# Plot Heston model asset prices and variance process
time <- seq(0, T, length.out = N + 1)

# Prepare data for plotting
df_asset <- data.frame(
  Time = rep(time, M),
  Price = as.vector(S_p),
  Simulation = rep(1:M, each = N + 1)
)

df_variance <- data.frame(
  Time = rep(time, M),
  Variance = as.vector(v_p),
  Simulation = rep(1:M, each = N + 1)
)

# Asset prices plot
p1 <- ggplot(df_asset, aes(x = Time, y = Price, group = Simulation)) +
  geom_line(alpha = 0.3) +
  labs(title = "Heston Model Asset Prices",
       x = "Time",
       y = "Asset Prices") +
  theme_minimal()

# Variance process plot
p2 <- ggplot(df_variance, aes(x = Time, y = Variance, group = Simulation)) +
  geom_line(alpha = 0.3) +
  labs(title = "Heston Model Variance Process",
       x = "Time",
       y = "Variance") +
  theme_minimal()

# Display plots
gridExtra::grid.arrange(p1, p2, ncol = 2)

# Simulate GBM process at time T
gbm <- S0 * exp((r - theta^2/2) * T + sqrt(theta) * sqrt(T) * rnorm(M))

# Prepare density data for comparison
df_density <- data.frame(
  Price = c(S_p[N + 1, ], S_n[N + 1, ], gbm),
  Model = factor(c(rep("rho = 0.98", M), 
                   rep("rho = -0.98", M), 
                   rep("GBM", M)))
)

# Density plot
ggplot(df_density, aes(x = Price, fill = Model, color = Model)) +
  geom_density(alpha = 0.5) +
  labs(title = "Asset Price Density under Heston Model",
       x = "S_T",
       y = "Density") +
  xlim(20, 180) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")

# Simulate with rho = -0.7
rho <- -0.7
result <- heston_model_sim(S0, v0, rho, kappa, theta, sigma, T, N, M)
S <- result$S

# Set strikes and complete MC option price for different strikes
K <- seq(20, 180, by = 2)

# Calculate put and call prices
puts <- sapply(K, function(k) exp(-r * T) * mean(pmax(k - S[N + 1, ], 0)))
calls <- sapply(K, function(k) exp(-r * T) * mean(pmax(S[N + 1, ] - k, 0)))

# Function to calculate Black-Scholes implied volatility
# Using a simple bisection method
black_scholes_price <- function(S, K, T, r, sigma, type = "call") {
  d1 <- (log(S/K) + (r + sigma^2/2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  if (type == "call") {
    price <- S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  } else if (type == "put") {
    price <- K * exp(-r * T) * pnorm(-d2) - S * pnorm(-d1)
  }
  return(price)
}

implied_volatility <- function(market_price, S, K, T, r, type, 
                               lower = 0.001, upper = 5, tol = 1e-6) {
  # Bisection method for implied volatility
  f <- function(sigma) {
    black_scholes_price(S, K, T, r, sigma, type) - market_price
  }
  
  if (f(lower) * f(upper) > 0) {
    return(NA)  # No solution found
  }
  
  for (i in 1:100) {
    mid <- (lower + upper) / 2
    f_mid <- f(mid)
    
    if (abs(f_mid) < tol) {
      return(mid)
    }
    
    if (f(lower) * f_mid < 0) {
      upper <- mid
    } else {
      lower <- mid
    }
  }
  return(mid)
}

# Calculate implied volatilities
put_ivs <- sapply(1:length(K), function(i) {
  implied_volatility(puts[i], S0, K[i], T, r, type = "put")
})

call_ivs <- sapply(1:length(K), function(i) {
  implied_volatility(calls[i], S0, K[i], T, r, type = "call")
})

# Create dataframe for plotting
df_iv <- data.frame(
  Strike = rep(K, 2),
  IV = c(call_ivs, put_ivs),
  Type = factor(rep(c("IV calls", "IV puts"), each = length(K)))
)

# Implied volatility smile plot
ggplot(df_iv, aes(x = Strike, y = IV, color = Type)) +
  geom_line(size = 1) +
  labs(title = "Implied Volatility Smile from Heston Model",
       x = "Strike",
       y = "Implied Volatility") +
  theme_minimal() +
  scale_color_manual(values = c("IV calls" = "blue", "IV puts" = "red"))

# Optional: Additional analysis
# Print summary statistics
cat("Summary Statistics for Terminal Asset Prices:\n")
cat("rho = 0.98: mean =", mean(S_p[N + 1, ]), "sd =", sd(S_p[N + 1, ]), "\n")
cat("rho = -0.98: mean =", mean(S_n[N + 1, ]), "sd =", sd(S_n[N + 1, ]), "\n")
cat("GBM: mean =", mean(gbm), "sd =", sd(gbm), "\n")