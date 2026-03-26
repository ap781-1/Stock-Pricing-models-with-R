# Enhanced Monte Carlo Option Pricing Model
# Features:
# 1. Multiple asset price models (GBM, Heston, Merton Jump Diffusion)
# 2. Control variates for variance reduction
# 3. Antithetic variates
# 4. Multiple exercise types (European, American via LSMC)
# 5. Confidence intervals
# 6. Convergence diagnostics
# 7. Error handling and validation

library(ggplot2)
library(dplyr)
library(tidyr)

# Main pricing function with multiple model options
enhanced_option_pricing <- function(nSim = 100000, 
                                    tau = 0.5, 
                                    r = 0.02, 
                                    sigma = 0.2, 
                                    S0 = 102, 
                                    K = 100,
                                    option_type = c("call", "put"),
                                    model = c("GBM", "Heston", "Merton"),
                                    exercise_type = c("European", "American"),
                                    use_control_variate = TRUE,
                                    use_antithetic = TRUE,
                                    n_steps = 252,  # for American options
                                    heston_params = list(kappa = 2.0, 
                                                        theta = 0.04, 
                                                        xi = 0.3, 
                                                        rho = -0.7,
                                                        v0 = 0.04),
                                    merton_params = list(lambda_j = 0.5,
                                                        mu_j = -0.05,
                                                        sigma_j = 0.3),
                                    confidence_level = 0.95,
                                    seed = NULL) {
  
  # Input validation
  if (!is.null(seed)) set.seed(seed)
  
  option_type <- match.arg(option_type)
  model <- match.arg(model)
  exercise_type <- match.arg(exercise_type)
  
  # Validate inputs
  if (any(c(nSim, tau, r, sigma, S0, K) <= 0)) {
    stop("All numerical parameters must be positive")
  }
  
  if (exercise_type == "American" && model != "GBM") {
    warning("American options only implemented for GBM model. Switching to GBM.")
    model <- "GBM"
  }
  
  # Price based on selected model
  if (exercise_type == "European") {
    result <- price_european(nSim, tau, r, sigma, S0, K, option_type, 
                            model, use_control_variate, use_antithetic,
                            heston_params, merton_params)
  } else {
    result <- price_american_lsm(nSim, tau, r, sigma, S0, K, option_type,
                                n_steps, use_antithetic, confidence_level)
  }
  
  # Add confidence intervals
  z_val <- qnorm((1 + confidence_level)/2)
  result$ci_lower <- result$price - z_val * result$sterr
  result$ci_upper <- result$price + z_val * result$sterr
  
  # Add relative error and efficiency metrics
  result$rel_error <- result$sterr / result$price
  result$efficiency <- 1 / (result$sterr^2 * result$nSim_effective)
  
  return(result)
}

# European option pricing with various models
price_european <- function(nSim, tau, r, sigma, S0, K, option_type,
                          model, use_control_variate, use_antithetic,
                          heston_params, merton_params) {
  
  if (use_antithetic) {
    nSim_half <- ceiling(nSim/2)
    Z <- rnorm(nSim_half)
    Z_anti <- -Z
    Z_all <- c(Z, Z_anti)
    nSim_effective <- length(Z_all)
  } else {
    Z_all <- rnorm(nSim)
    nSim_effective <- nSim
  }
  
  # Generate terminal prices based on selected model
  if (model == "GBM") {
    ST <- generate_gbm_paths(S0, r, sigma, tau, Z_all)
  } else if (model == "Heston") {
    ST <- generate_heston_paths(S0, r, tau, Z_all, heston_params, nSim_effective)
  } else if (model == "Merton") {
    ST <- generate_merton_paths(S0, r, sigma, tau, Z_all, merton_params)
  }
  
  # Calculate payoffs
  if (option_type == "call") {
    payoffs_raw <- pmax(ST - K, 0)
  } else {
    payoffs_raw <- pmax(K - ST, 0)
  }
  
  payoffs <- exp(-r * tau) * payoffs_raw
  
  # Control variate using Black-Scholes price
  if (use_control_variate && model == "GBM") {
    bs_price <- black_scholes(S0, K, tau, r, sigma, option_type)
    
    # Use discounted asset price as control variate
    control_variate <- exp(-r * tau) * ST
    
    # Calculate optimal coefficient
    cov_xy <- cov(payoffs, control_variate)
    var_x <- var(control_variate)
    beta <- cov_xy / var_x
    
    # Adjust payoffs
    control_mean <- S0  # E[discounted asset price] = S0
    payoffs_adjusted <- payoffs - beta * (control_variate - control_mean)
    
    price <- mean(payoffs_adjusted)
    sterr <- sd(payoffs_adjusted) / sqrt(nSim_effective)
    
    # Store beta for reporting
    control_beta <- beta
  } else {
    price <- mean(payoffs)
    sterr <- sd(payoffs) / sqrt(nSim_effective)
    control_beta <- NULL
  }
  
  return(list(
    price = price,
    sterr = sterr,
    nSim_effective = nSim_effective,
    model = model,
    option_type = option_type,
    control_beta = control_beta,
    raw_payoffs = payoffs
  ))
}

# Geometric Brownian Motion paths
generate_gbm_paths <- function(S0, r, sigma, tau, Z) {
  drift <- (r - 0.5 * sigma^2) * tau
  diffusion <- sigma * sqrt(tau) * Z
  return(S0 * exp(drift + diffusion))
}

# Heston model paths (simplified - Euler discretization)
generate_heston_paths <- function(S0, r, tau, Z, params, nSim) {
  kappa <- params$kappa
  theta <- params$theta
  xi <- params$xi
  rho <- params$rho
  v0 <- params$v0
  
  n_steps <- 252  # Daily steps for Heston
  dt <- tau / n_steps
  
  # Generate correlated random numbers
  Z1 <- matrix(rnorm(nSim * n_steps), nrow = nSim)
  Z2 <- rho * Z1 + sqrt(1 - rho^2) * matrix(rnorm(nSim * n_steps), nrow = nSim)
  
  S <- matrix(S0, nrow = nSim, ncol = 1)
  v <- matrix(v0, nrow = nSim, ncol = 1)
  
  for (i in 1:n_steps) {
    # Euler discretization for variance (ensure non-negative)
    v_next <- v[, i] + kappa * (theta - v[, i]) * dt + 
              xi * sqrt(v[, i] * dt) * Z2[, i]
    v_next <- pmax(v_next, 0)  # Absorption at zero
    
    # Euler for asset price
    S_next <- S[, i] * exp((r - 0.5 * v[, i]) * dt + 
                           sqrt(v[, i] * dt) * Z1[, i])
    
    S <- cbind(S, S_next)
    v <- cbind(v, v_next)
  }
  
  return(S[, n_steps + 1])
}

# Merton Jump Diffusion model
generate_merton_paths <- function(S0, r, sigma, tau, Z, params) {
  lambda <- params$lambda_j
  mu_j <- params$mu_j
  sigma_j <- params$sigma_j
  
  # Number of jumps (Poisson)
  n_jumps <- rpois(length(Z), lambda * tau)
  
  # Jump sizes (log-normal)
  jump_sizes <- sapply(n_jumps, function(n) {
    if (n > 0) {
      sum(rnorm(n, mu_j, sigma_j))
    } else {
      0
    }
  })
  
  # Diffusion component
  diffusion <- (r - 0.5 * sigma^2 - lambda * (exp(mu_j + 0.5 * sigma_j^2) - 1)) * tau +
               sigma * sqrt(tau) * Z
  
  return(S0 * exp(diffusion + jump_sizes))
}

# Black-Scholes formula for control variate
black_scholes <- function(S0, K, T, r, sigma, option_type) {
  d1 <- (log(S0/K) + (r + sigma^2/2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  if (option_type == "call") {
    price <- S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  } else {
    price <- K * exp(-r * T) * pnorm(-d2) - S0 * pnorm(-d1)
  }
  
  return(price)
}

# American option pricing using Longstaff-Schwartz (LSM) method
price_american_lsm <- function(nSim, tau, r, sigma, S0, K, option_type,
                              n_steps, use_antithetic, confidence_level) {
  
  dt <- tau / n_steps
  discount <- exp(-r * dt)
  
  # Generate paths with antithetic variates
  if (use_antithetic) {
    nSim_half <- ceiling(nSim/2)
    Z <- matrix(rnorm(nSim_half * n_steps), nrow = nSim_half, ncol = n_steps)
    Z_anti <- -Z
    Z_all <- rbind(Z, Z_anti)
    nSim_effective <- nrow(Z_all)
  } else {
    Z_all <- matrix(rnorm(nSim * n_steps), nrow = nSim, ncol = n_steps)
    nSim_effective <- nSim
  }
  
  # Generate price paths
  S <- matrix(0, nrow = nSim_effective, ncol = n_steps + 1)
  S[, 1] <- S0
  
  for (i in 1:n_steps) {
    S[, i + 1] <- S[, i] * exp((r - 0.5 * sigma^2) * dt + 
                                sigma * sqrt(dt) * Z_all[, i])
  }
  
  # Initialize cash flow matrix
  if (option_type == "call") {
    payoff <- pmax(S - K, 0)
  } else {
    payoff <- pmax(K - S, 0)
  }
  
  cash_flow <- payoff[, n_steps + 1]
  
  # Backward induction
  for (i in n_steps:1) {
    in_the_money <- payoff[, i] > 0
    
    if (sum(in_the_money) > 3) {  # Need enough points for regression
      # Regression of discounted future cash flows on current price
      X <- S[in_the_money, i]
      Y <- cash_flow[in_the_money] * discount^(n_steps - i + 1)
      
      # Basis functions (polynomials)
      basis <- cbind(1, X, X^2, X^3)
      reg <- lm(Y ~ X + I(X^2) + I(X^3) - 1, data = data.frame(X = X))
      
      # Expected continuation value
      continuation <- predict(reg, newdata = data.frame(X = X))
      
      # Exercise decision
      exercise <- payoff[in_the_money, i] > continuation
      
      # Update cash flows
      cash_flow[in_the_money][exercise] <- payoff[in_the_money, i][exercise]
    }
  }
  
  # Discount to present value
  option_price <- mean(cash_flow * discount^(n_steps))
  option_sterr <- sd(cash_flow * discount^(n_steps)) / sqrt(nSim_effective)
  
  return(list(
    price = option_price,
    sterr = option_sterr,
    nSim_effective = nSim_effective,
    model = "LSM (American)",
    option_type = option_type
  ))
}

# Convergence analysis function
convergence_analysis <- function(nSim_values, tau, r, sigma, S0, K, 
                                 option_type = "call", n_reps = 10) {
  
  results <- data.frame()
  
  for (n in nSim_values) {
    prices <- numeric(n_reps)
    sterrs <- numeric(n_reps)
    
    for (i in 1:n_reps) {
      res <- enhanced_option_pricing(nSim = n, tau = tau, r = r, 
                                     sigma = sigma, S0 = S0, K = K,
                                     option_type = option_type,
                                     use_control_variate = TRUE,
                                     use_antithetic = TRUE)
      prices[i] <- res$price
      sterrs[i] <- res$sterr
    }
    
    results <- rbind(results, data.frame(
      nSim = n,
      mean_price = mean(prices),
      sd_price = sd(prices),
      mean_sterr = mean(sterrs),
      rel_error = mean(sterrs) / mean(prices)
    ))
  }
  
  return(results)
}

# Visualization function
plot_convergence <- function(convergence_results) {
  p1 <- ggplot(convergence_results, aes(x = nSim, y = mean_price)) +
    geom_line() +
    geom_ribbon(aes(ymin = mean_price - sd_price, 
                    ymax = mean_price + sd_price), alpha = 0.3) +
    scale_x_log10() +
    labs(title = "Convergence of Monte Carlo Option Price",
         x = "Number of Simulations (log scale)",
         y = "Option Price") +
    theme_minimal()
  
  p2 <- ggplot(convergence_results, aes(x = nSim, y = rel_error)) +
    geom_line(color = "red") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Relative Error Convergence",
         x = "Number of Simulations (log scale)",
         y = "Relative Error (log scale)") +
    theme_minimal()
  
  return(list(price_plot = p1, error_plot = p2))
}

# Example usage
set.seed(123)

# European call with GBM and control variates
result1 <- enhanced_option_pricing(
  nSim = 100000,
  tau = 0.5,
  r = 0.02,
  sigma = 0.2,
  S0 = 102,
  K = 100,
  option_type = "call",
  model = "GBM",
  use_control_variate = TRUE,
  use_antithetic = TRUE
)

print(result1)

# European put with Heston model
result2 <- enhanced_option_pricing(
  nSim = 50000,
  tau = 0.5,
  r = 0.02,
  sigma = 0.2,
  S0 = 102,
  K = 100,
  option_type = "put",
  model = "Heston",
  use_control_variate = FALSE
)

print(result2)

# American put using Longstaff-Schwartz
result3 <- enhanced_option_pricing(
  nSim = 20000,
  tau = 0.5,
  r = 0.02,
  sigma = 0.2,
  S0 = 102,
  K = 100,
  option_type = "put",
  exercise_type = "American",
  n_steps = 50
)

print(result3)

# Convergence analysis
nSim_values <- c(1000, 5000, 10000, 50000, 100000)
conv_results <- convergence_analysis(nSim_values, tau = 0.5, r = 0.02, 
                                     sigma = 0.2, S0 = 102, K = 100)
print(conv_results)

# Plot convergence
plots <- plot_convergence(conv_results)
print(plots$price_plot)
print(plots$error_plot)
