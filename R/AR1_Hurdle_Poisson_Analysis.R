# --- Packages ----
#if(!require("pacman")) install.packages("pacman")

pacman::p_load(zoo, xts, forecast, tseries, ggplot2, pscl, stats4, numDeriv, xtable)

# --- 1) Read and process ----
df <- read.csv("Data/WHO_Kenya_daily_deaths.csv", stringsAsFactors = FALSE)
df$date <- as.Date(df$date, format="%d/%m/%Y")
df <- df[order(df$date), ]
df$deaths <- as.integer(df$deaths)
rownames(df) <- NULL

# create lagged variable (must be done before any summaries using it)
df$deaths_lag1 <- c(NA, head(df$deaths, -1))

# Quick checks
cat("Rows:", nrow(df), " NA in deaths:", sum(is.na(df$deaths)), " NA in lag:", sum(is.na(df$deaths_lag1)), "\n")


######## Summary statistics, dispersion and zero proportion ########

min_deaths  <- min(df$deaths)
q1_deaths   <- quantile(df$deaths, 0.25)
median_deaths <- median(df$deaths)
mean_deaths <- mean(df$deaths)
q3_deaths   <- quantile(df$deaths, 0.75)
max_deaths  <- max(df$deaths)

var_deaths  <- var(df$deaths)
dispersion_index <- var_deaths / mean_deaths
prop_zero   <- mean(df$deaths == 0)

# Lag variable

min_lag  <- min(df$deaths_lag1, na.rm=TRUE)
q1_lag   <- quantile(df$deaths_lag1, 0.25, na.rm=TRUE)
median_lag <- median(df$deaths_lag1, na.rm=TRUE)
mean_lag <- mean(df$deaths_lag1, na.rm=TRUE)
q3_lag   <- quantile(df$deaths_lag1, 0.75, na.rm=TRUE)
max_lag  <- max(df$deaths_lag1, na.rm=TRUE)

# Print diagnostics
mean_deaths
var_deaths
dispersion_index
prop_zero



summary_table <- data.frame(
  Statistic = c("Minimum",
                "1st Quartile",
                "Median",
                "Mean",
                "3rd Quartile",
                "Maximum",
                "Variance",
                "Dispersion Index (Var/Mean)",
                "Proportion of zeros"),
  
  "Deaths $(X_t)$" = c(
    min(df$deaths),
    as.numeric(quantile(df$deaths,0.25)),
    median(df$deaths),
    mean(df$deaths),
    as.numeric(quantile(df$deaths,0.75)),
    max(df$deaths),
    var(df$deaths),
    var(df$deaths)/mean(df$deaths),
    mean(df$deaths==0)
  ),
  
  "Lagged Deaths $(X_{t-1})$" = c(
    min(df$deaths_lag1, na.rm=TRUE),
    as.numeric(quantile(df$deaths_lag1,0.25, na.rm=TRUE)),
    median(df$deaths_lag1, na.rm=TRUE),
    mean(df$deaths_lag1, na.rm=TRUE),
    as.numeric(quantile(df$deaths_lag1,0.75, na.rm=TRUE)),
    max(df$deaths_lag1, na.rm=TRUE),
    NA,
    NA,
    NA
  ),
  
  check.names = FALSE
)

print(
  xtable(summary_table,
         caption = "Summary statistics of daily COVID-19 deaths in Kenya",
         label = "tab:summary_stats",
         digits = 3),
  include.rownames = FALSE,
  sanitize.text.function = identity,
  caption.placement = "top",
  file = "summary_stats.tex"
)


# plot daily deaths
ggplot(df, aes(x=date, y=deaths)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Deaths",
    title = "Daily Deaths in Kenya"
  )

# Save as PDF

ggsave("daily_deaths.pdf", width=7, height=5)

# OR save as PNG
#ggsave("daily_deaths.png", width=7, height=5, dpi=300)


ggplot(df, aes(x=date)) +
  geom_line(aes(y=deaths, color="Deaths")) +
  geom_line(aes(y=deaths_lag1, color="Lagged Deaths (1-day)")) +
  labs(x="Date", y="Deaths", title="Daily deaths and 1-day lagged deaths") +
  scale_color_manual(values=c("Deaths"="blue", "Lagged Deaths (1-day)"="red"))

ggsave("deaths_with_lag.pdf", width=7, height=5)

# --- ) EDA: histogram, ACF/PACF ----
hist(df$deaths, breaks=40, main="Histogram of daily deaths", xlab="Deaths")
acf(df$deaths, na.action = na.contiguous, main="ACF (counts)")
pacf(df$deaths, na.action = na.contiguous, main="PACF (counts)")

# Show in RStudio
#par(mar=c(5,5,4,2))
#hist(df$deaths, breaks=40, main="Histogram of daily deaths", xlab="Deaths")

# Save histogram
png("histogram_daily_deaths.png", width=7, height=5, units="in", res=300)
par(mar=c(5,5,4,2))  
hist(df$deaths, breaks=40, main="Histogram of daily deaths", xlab="Deaths")
dev.off()

# Save ACF
png("acf_daily_deaths.png", width=7, height=5, units="in", res=300)
acf(df$deaths, na.action=na.contiguous, main="ACF (counts)")
dev.off()

# Save PACF
png("pacf_daily_deaths.png", width=7, height=5, units="in", res=300)
pacf(df$deaths, na.action=na.contiguous, main="PACF (counts)")
dev.off()

# zero/positive separation
df$zero <- ifelse(df$deaths == 0, 1, 0)
prop_zero <- mean(df$zero, na.rm=TRUE)
cat("Proportion zeros:", round(prop_zero,3), "\n")
acf(df$zero, main="ACF (zero indicator)")

# Save ACF of zero-indicator
png("acf_zero.png", width=800, height=600)
acf(df$zero, main="ACF (zero indicator)")
dev.off()


# PACF for zero-indicator
pacf(df$zero, main = "PACF (zero indicator)")

# Save PACF plot for zero-indicator
png("pacf_zero.png", width=800, height=600)
pacf(df$zero, main = "PACF (zero indicator)")
dev.off()

########## Take note- New revised Zero Indicator Plots

# zero indicator
df$zero <- ifelse(df$deaths == 0, 1, 0)

# ------------------ ACF ------------------
acf_data <- Acf(df$zero, lag.max = 30, plot = FALSE)

acf_df <- with(acf_data, data.frame(lag, acf))

p1 <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "orange", linewidth = 1) +
  geom_segment(aes(xend = lag, yend = 0), color = "orange", linewidth = 1.2) +
  geom_point(color = "orange", size = 3) +
  geom_ribbon(aes(ymin = -1.96/sqrt(length(df$zero)),
                  ymax = 1.96/sqrt(length(df$zero))),
              fill = "#F5C97A", alpha = 0.5) +
  labs(title = "Autocorrelation Function (ACF) of Zero Indicator",
       x = "Lag",
       y = "ACF") +
  theme_minimal(base_size = 16)

ggsave("acf_zero.png", p1, width = 10, height = 6)


# ------------------ PACF ------------------
pacf_data <- Pacf(df$zero, lag.max = 30, plot = FALSE)

pacf_df <- with(pacf_data, data.frame(lag = lag, pacf = acf))

p2 <- ggplot(pacf_df, aes(x = lag, y = pacf)) +
  geom_hline(yintercept = 0, color = "orange", linewidth = 1) +
  geom_segment(aes(xend = lag, yend = 0), color = "orange", linewidth = 1.2) +
  geom_point(color = "orange", size = 3) +
  geom_ribbon(aes(ymin = -1.96/sqrt(length(df$zero)),
                  ymax = 1.96/sqrt(length(df$zero))),
              fill = "#F5C97A", alpha = 0.5) +
  labs(title = "Partial Autocorrelation Function (PACF) of Zero Indicator",
       x = "Lag",
       y = "PACF") +
  theme_minimal(base_size = 16)

ggsave("pacf_zero.png", p2, width = 10, height = 6)

#############################################################################
library(dplyr)
library(ggplot2)
# --- Prepare modelling dataset ----

# For modelling we will use only rows where 'deaths' is not NA and lag is available
df_model <- df %>% select(date, deaths, deaths_lag1) %>% na.omit()
rownames(df_model) <- NULL

cat("Rows used for modeling (complete cases):", nrow(df_model), "\n")

# --- 2) Starting values via simplified regressions (use complete cases) ----
# Zero part (logit) using indicator of zero
df_model$zero <- ifelse(df_model$deaths == 0, 1, 0)
glm_zero <- glm(zero ~ deaths_lag1, family = binomial(link = "logit"), data = df_model)
start_alpha <- as.numeric(coef(glm_zero)[1])
start_eta   <- ifelse(length(coef(glm_zero)) >= 2, as.numeric(coef(glm_zero)[2]), 0)

# Positive part (Poisson) on positive counts only
df_pos <- df_model %>% filter(deaths > 0)
if(nrow(df_pos) >= 5){
  glm_pos <- glm(deaths ~ deaths_lag1, family = poisson(link = "log"), data = df_pos)
  start_gamma <- as.numeric(coef(glm_pos)[1])
  start_theta <- ifelse(length(coef(glm_pos)) >= 2, as.numeric(coef(glm_pos)[2]), 0)
} else {
  # fallback small defaults if too few positives
  start_gamma <- log(max(1, mean(df_model$deaths[df_model$deaths>0], na.rm=TRUE)))
  start_theta <- 0.05
}

# psi is a slope in the log-linear mean (treat as real-valued)
start_psi <- 0.1

# Bundle starting values: NOTE -- psi is on the real line (no logit)
start_par <- c(
  alpha = start_alpha,
  eta   = start_eta,
  beta  = 0.1,           # small default for direct lag on zero-prob
  gamma = start_gamma,
  theta = start_theta,
  psi   = start_psi
)
print(start_par)

# --- 3) Stabilized negative log-likelihood (robust) ----
neg_logLik_AR1_HP <- function(par, data){
  # Accept named or unnamed 'par'
  if(is.null(names(par))){
    # assume order: alpha, eta, beta, gamma, theta, psi
    names(par) <- c("alpha","eta","beta","gamma","theta","psi")
  }
  # unpack (psi is plain slope, NOT a probability)
  alpha <- as.numeric(par["alpha"])
  eta   <- as.numeric(par["eta"])
  beta  <- as.numeric(par["beta"])
  gamma <- as.numeric(par["gamma"])
  theta <- as.numeric(par["theta"])
  psi   <- as.numeric(par["psi"])   # <- NO plogis()
  
  # Data vectors (assume data is complete cases with lag available)
  deaths <- as.numeric(data$deaths)
  n <- length(deaths)
  if(n <= 1) return(1e10)
  # Initialization
  pi_prev <- mean(deaths == 0, na.rm = TRUE)
  mu_prev <- mean(deaths[deaths > 0], na.rm = TRUE)
  if(is.na(mu_prev) || mu_prev <= 0) mu_prev <- 1.0
  
  loglik <- 0.0
  
  for(t in 1:n){
    x_t <- deaths[t]
    # if missing value encountered, skip (shouldn't happen because we used na.omit above),
    # but protect anyway
    if(is.na(x_t)) next
    
    x_lag <- if(t==1) 0 else deaths[t-1]
    
    # Zero probability (logit link). Use continuous feedback via pi_prev (not boolean).
    linear_pi <- alpha + eta * pi_prev + beta * x_lag
    # clamp to avoid overflow
    linear_pi <- pmin(pmax(linear_pi, -30), 30)
    pi_t <- 1 / (1 + exp(-linear_pi))
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8)
    
    # Poisson mean (log link)
    linear_mu <- gamma + theta * log(mu_prev + 1e-8) + psi * x_lag   # use log(mu_prev) to stabilize scale
    linear_mu <- pmin(pmax(linear_mu, -30), 30)
    mu_t <- exp(linear_mu)
    mu_t <- pmax(mu_t, 1e-8)
    
    # truncated Poisson denominator for positive part (1 - P(Y=0))
    denom <- 1 - exp(-mu_t)
    denom <- pmax(denom, 1e-12)
    
    # likelihood for this t
    if(x_t == 0){
      p <- pi_t
    } else {
      poisson_p_log <- dpois(x_t, lambda = mu_t, log = TRUE)
      # compute log p then exp if necessary: but we aggregate log-lik
      log_p <- log(1 - pi_t) + poisson_p_log - log(denom)
      p <- exp(log_p)
      # protect
      if(!is.finite(p) || p <= 0) p <- 1e-12
    }
    # accumulate log-likelihood
    loglik <- loglik + log(p)
    
    # update recursive states
    pi_prev <- pi_t
    mu_prev <- mu_t
  }
  
  # return negative log-likelihood (optimizer minimizes)
  neg_ll <- -loglik
  # protect: ensure finite numeric
  if(!is.finite(neg_ll) || is.na(neg_ll)) return(1e10)
  return(neg_ll)
}

# --- 4) Test the log-likelihood at starting values ----
test_val <- tryCatch({
  neg_logLik_AR1_HP(start_par, df_model)
}, error = function(e){
  cat("Error evaluating neg-logLik at start_par:", conditionMessage(e), "\n")
  return(1e10)
})
cat("neg-logLik at start:", test_val, "\n")

# --- CHECK POINTS) Benchmark thesis estimates against stabilized likelihood ----

thesis_par <- c(
  alpha = -1.70866,
  eta   =  3.98443,
  beta  = -0.26441,
  gamma =  0.50976,
  theta =  0.62702,
  psi   =  0.02810
)

thesis_nll <- neg_logLik_AR1_HP(
  thesis_par,
  data = df_model
)

cat(
  "NLL at thesis estimates:",
  thesis_nll,
  "\n"
)

cat(
  "Is the value finite?",
  is.finite(thesis_nll),
  "\n"
)

# --- CHECK POINTS) Optimization starting from thesis estimates ----

opt_thesis_start <- optim(
  par = thesis_par,
  fn = function(par) neg_logLik_AR1_HP(par, data = df_model),
  method = "BFGS",
  control = list(maxit = 5000, trace = 1)
)

print(opt_thesis_start)

cat("\nStarting NLL:", thesis_nll, "\n")
cat("Optimized NLL:", opt_thesis_start$value, "\n")
cat("Convergence code:", opt_thesis_start$convergence, "\n")

cat("\nOptimized parameters:\n")
print(opt_thesis_start$par)

# --- CHECK POINTS) Hessian, Standard Errors, z-values and p-values ----

# Use the confirmed thesis solution
est <- opt_thesis_start$par

# Compute numerical Hessian of the stabilized negative log-likelihood
hess <- numDeriv::hessian(
  func = function(par) {
    neg_logLik_AR1_HP(par, data = df_model)
  },
  x = est
)

# Variance-covariance matrix
vcov_est <- tryCatch(
  solve(hess),
  error = function(e) {
    stop("Hessian could not be inverted: ", conditionMessage(e))
  }
)

# Standard errors
se <- sqrt(diag(vcov_est))

# z-values
z_values <- est / se
# Two-sided p-values
p_values <- 2 * (1 - pnorm(abs(z_values)))

# Results table
res_table <- data.frame(
  Parameter = names(est),
  Estimate  = as.numeric(est),
  StdErr    = as.numeric(se),
  z_value   = as.numeric(z_values),
  p_value   = as.numeric(p_values),
  row.names = NULL
)

# Display with sufficient precision for comparison
print(res_table)



# --- 5) Optimization: try nlminb first, then optim fallbacks ----
# Parameters are estimated on the real line.
# No parameter bounds are imposed in the primary nlminb optimization.

lower_bounds <- rep(-Inf, length(start_par))
upper_bounds <- rep(Inf, length(start_par))
names(lower_bounds) <- names(start_par); names(upper_bounds) <- names(start_par)

# Primary: nlminb (often more stable for unconstrained real params)
opt_try <- tryCatch({
  nlminb(start = start_par,
         objective = function(p) neg_logLik_AR1_HP(p, df_model),
         control = list(iter.max = 20000, eval.max = 20000))
}, error = function(e){
  warning("nlminb failed: ", conditionMessage(e))
  return(NULL)
})

if(is.null(opt_try) || !is.numeric(opt_try$objective)){
  # fallback to optim (BFGS)
  cat("nlminb failed or returned invalid objective; trying optim (BFGS)...\n")
  opt_try2 <- tryCatch({
    optim(par = start_par,
          fn = function(p) neg_logLik_AR1_HP(p, df_model),
          method = "BFGS",
          control = list(maxit = 20000))
  }, error = function(e){
    warning("optim (BFGS) failed: ", conditionMessage(e))
    return(NULL)
  })
  opt_final <- opt_try2
} else {
  opt_final <- opt_try
}

# If optim produced 'par' in $par or $estimate, normalize to get 'est'
if(is.null(opt_final)){
  stop("Both optimizers failed. Inspect input data and starting values.")
}
est <- if(!is.null(opt_final$par)) opt_final$par else if(!is.null(opt_final$estimate)) opt_final$estimate else stop("Cannot extract estimates")

cat("Convergence info (nlminb/optim object):\n"); print(opt_final)

# --- 6a) Hessian (numDeriv) and likelihood-based variance–covariance matrix ---- ----
hess_use <- tryCatch({
  numDeriv::hessian(func = function(p) neg_logLik_AR1_HP(p, df_model), x = est)
}, error = function(e){
  warning("Hessian computation failed: ", conditionMessage(e))
  return(matrix(NA, nrow=length(est), ncol=length(est)))
})

# Try to invert Hessian to get variance-covariance
vcov_est <- tryCatch({
  solve(hess_use)
}, error = function(e){
  warning("Hessian inversion failed: ", conditionMessage(e))
  return(matrix(NA, nrow=length(est), ncol=length(est)))
})

# Standard errors (with names for alignment)
se <- rep(NA, length(est))
if(all(!is.na(diag(vcov_est)))){
  se <- sqrt(pmax(0, diag(vcov_est)))
}
names(se) <- names(est)

# --- Make sure parameter names align and build results table (no special logit_psi handling) ---
# 'est' is a named numeric vector; 'se' likewise

# Guard against division by zero/NA
z_vals <- rep(NA, length(est))
p_vals <- rep(NA, length(est))
names(z_vals) <- names(est)
names(p_vals) <- names(est)

valid_se <- !is.na(se) & se > 0
z_vals[valid_se] <- est[valid_se] / se[valid_se]
p_vals[valid_se] <- 2 * (1 - pnorm(abs(z_vals[valid_se])))

res_table <- data.frame(
  Parameter = names(est),
  Estimate  = round(as.numeric(est), 5),
  StdErr    = round(as.numeric(se), 5),
  z_value   = round(z_vals, 3),
  p_value = ifelse(p_vals < 0.001, "<0.001", signif(p_vals, 3)),
  row.names = NULL,
  stringsAsFactors = FALSE
)

print(res_table)


# --- Export to LaTeX with xtable ---
xt <- xtable(res_table,
             digits = c(0, 0, 5, 5, 3, 3),
             caption = "Parameter Estimates (AR(1) Hurdle-Poisson) with SEs, z and p-values",
             label   = "tab:par_estimates")

print(xt,
      include.rownames = FALSE,  # ✅ no duplicate column
      caption.placement = "top",
      file = "results_table.tex")

# --- 6b) Extract model fit statistics (Goodness-of-Fit) ---

# Compute log-likelihood at the estimated parameters
ll_val <- -neg_logLik_AR1_HP(est, df_model)   # since neg_logLik returns NEGATIVE log-lik

# Number of free parameters (length of estimated vector, including logit_psi)
k <- length(est)
n <- nrow(df_model)


# AIC and BIC
aic_val <- -2 * ll_val + 2 * k
bic_val <- -2 * ll_val + log(n) * k

cat("Log-likelihood:", ll_val, "\n")
cat("AIC:", aic_val, "\n")
cat("BIC:", bic_val, "\n")

# Create small data.frame for LaTeX export
fit_stats <- data.frame(
  Criterion = c("Log-Likelihood", "AIC", "BIC"),
  Value     = c(round(ll_val, 2),
                round(aic_val, 2),
                round(bic_val, 2))
)

# Export to LaTeX
xt_fit <- xtable(fit_stats,
                 caption = "Model Fit Statistics for AR(1) Hurdle-Poisson model",
                 label   = "tab:fit_stats")

print(xt_fit,
      include.rownames = FALSE,
      caption.placement = "top",
      file = "fit_stats.tex")

############### --- Baseline Hurdle Poisson (no AR terms)

neg_loglik_hp <- function(par, data){
  alpha <- par[1]
  gamma <- par[2]
  
  deaths <- data$deaths
  n <- length(deaths)
  
  loglik <- 0
  
  # Constant zero probability and positive-count parameter
  pi <- 1 / (1 + exp(-alpha))
  mu <- exp(gamma)
  
  for(t in 1:n){
    
    if(deaths[t] == 0){
      loglik <- loglik + log(pi + 1e-10)
    } else {
      loglik <- loglik +
        log(1 - pi + 1e-10) +
        dpois(deaths[t], mu, log=TRUE) -
        log(1 - dpois(0, mu) + 1e-10)
    }
  }
  
  return(-loglik)
}

# Initial values
init_hp <- c(alpha=-1, gamma=0)

# Fit model
fit_hp <- optim(
  init_hp,
  neg_loglik_hp,
  data=df_model,
  method="BFGS",
  hessian=TRUE
)




# Extract log-likelihood
loglik_hp <- -fit_hp$value

# Number of parameters
k_hp <- length(fit_hp$par)

# Sample size
n <- nrow(df_model)

# AIC and BIC
AIC_hp <- -2 * loglik_hp + 2 * k_hp
BIC_hp <- -2 * loglik_hp + log(n) * k_hp

######## --- Combine with AR(1) results

comparison_stats <- data.frame(
  Model = c("Hurdle Poisson", "AR(1) Hurdle Poisson"),
  LogLikelihood = c(round(loglik_hp, 2), round(ll_val, 2)),
  AIC = c(round(AIC_hp, 2), round(aic_val, 2)),
  BIC = c(round(BIC_hp, 2), round(bic_val, 2))
)

# --- EXPORT TO LATEX ---


xt_comp <- xtable(comparison_stats,
                  caption = "Comparison of baseline and AR(1) Hurdle Poisson models",
                  label   = "tab:comparison",
                  digits  = c(0, 0, 2, 2, 2))

print(xt_comp,
      include.rownames = FALSE,
      caption.placement = "top",
      file = "comparison_table.tex")




################################################################################

# ---7) Function to compute fitted pi_t (zero probability) and mu_t (positive mean)
# 7.1 Compute Fitted Probabilities and Expected Counts

compute_fitted_terms <- function(par, data){
  alpha <- par["alpha"]
  eta   <- par["eta"]
  beta  <- par["beta"]
  gamma <- par["gamma"]
  theta <- par["theta"]
  psi   <- par["psi"]   # plain coefficient
  
  deaths <- data$deaths
  n <- length(deaths)
  
  pi_vec <- numeric(n)
  mu_vec <- numeric(n)
  
  pi_prev <- mean(deaths == 0)
  mu_prev <- mean(deaths[deaths > 0])
  if(is.na(mu_prev) || mu_prev <= 0) mu_prev <- 1.0
  
  for(t in 1:n){
    x_lag <- if(t == 1) 0 else deaths[t-1]
    linear_pi <- as.numeric(alpha) + as.numeric(eta) * pi_prev + as.numeric(beta) * x_lag
    linear_pi <- pmin(pmax(linear_pi, -30), 30)
    pi_t <- 1 / (1 + exp(-linear_pi))
    pi_vec[t] <- pi_t
    
    linear_mu <- as.numeric(gamma) + as.numeric(theta) * log(mu_prev + 1e-8) + as.numeric(psi) * x_lag
    linear_mu <- pmin(pmax(linear_mu, -30), 30)
    mu_t <- exp(linear_mu)
    mu_vec[t] <- mu_t
    
    pi_prev <- pi_t
    mu_prev <- mu_t
  }
  
  return(list(pi_t = pi_vec, mu_t = mu_vec))
}

# Compute fitted terms
fitted_terms <- compute_fitted_terms(est, df_model)
df_model$pi_hat <- fitted_terms$pi_t
df_model$mu_hat <- fitted_terms$mu_t



xtable(df_model[1:10, c("date","deaths","pi_hat","mu_hat")])


#################################################################################
##### Revised Residual Analysis####


# --- 7.2 Residual Analysis ---

# Conditional mean and variance of the fitted Hurdle-Poisson distribution
p_zero <- df_model$pi_hat
mu <- df_model$mu_hat

# Probability of a positive count under the hurdle model
p_positive <- 1 - p_zero

# Zero-truncated Poisson denominator
trunc_prob <- -expm1(-mu)   # equivalent to 1 - exp(-mu), numerically stable

# Conditional mean of the full Hurdle-Poisson distribution
hurdle_mean <- p_positive * mu / trunc_prob

# Conditional second moment of the full Hurdle-Poisson distribution
hurdle_second_moment <- p_positive * mu * (1 + mu) / trunc_prob

# Conditional variance of the full Hurdle-Poisson distribution
hurdle_variance <- hurdle_second_moment - hurdle_mean^2

# Numerical protection
hurdle_variance <- pmax(hurdle_variance, 1e-12)

# Pearson residuals
df_model$pearson_res <- 
  (df_model$deaths - hurdle_mean) / sqrt(hurdle_variance)

# Plot residuals over time
library(ggplot2)

ggplot(df_model, aes(x = date, y = pearson_res)) +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Pearson Residuals Over Time",
    y = "Pearson Residuals",
    x = "Date"
  ) +
  theme_minimal()

ggsave("pearson_residuals.pdf")

# --- ACF of Pearson residuals ---

pdf("acf_pearson_residuals.pdf")

acf(
  df_model$pearson_res,
  main = "ACF of Pearson Residuals",
  lag.max = 30
)

dev.off()

#######CHECK POINTS FOR THE ABOVE

#1. Numerical ACF of the corrected Pearson residuals
acf_res <- acf(
  df_model$pearson_res,
  lag.max = 30,
  plot = FALSE
)

acf_values <- as.numeric(acf_res$acf)[-1]

acf_table <- data.frame(
  Lag = 1:30,
  ACF = round(acf_values, 4)
)

print(acf_table)

#2. Calculate the approximate 95% limits

n_res <- length(na.omit(df_model$pearson_res))

acf_limit <- 1.96 / sqrt(n_res)

cat("Number of residuals:", n_res, "\n")
cat("Approximate 95% ACF limits: +/-", round(acf_limit, 4), "\n")

#3. Get a basic residual summary

summary(df_model$pearson_res)

cat(
  "Mean Pearson residual:",
  mean(df_model$pearson_res, na.rm = TRUE),
  "\n"
)

cat(
  "SD Pearson residual:",
  sd(df_model$pearson_res, na.rm = TRUE),
  "\n"
)


#4.

Box.test(
df_model$pearson_res,
lag = 10,
type = "Ljung-Box"
)

#5.

Box.test(
df_model$pearson_res,
lag = 20,
type = "Ljung-Box"
) 

#6.

Box.test(
  df_model$pearson_res,
  lag = 30,
  type = "Ljung-Box"
)

############################################
# Old 7.2 Residual Analysis

# Pearson residuals
df_model$pearson_res <- (df_model$deaths - df_model$mu_hat) / sqrt(df_model$mu_hat * (1 - df_model$pi_hat))

# Plot residuals over time
library(ggplot2)
ggplot(df_model, aes(x=date, y=pearson_res)) +
  geom_line(color="steelblue") +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  labs(title="Pearson Residuals Over Time", y="Pearson Residuals", x="Date") +
  theme_minimal()

ggsave("pearson_residuals.pdf")

# --- Compute ACF of Pearson residuals and save ---
pdf("acf_pearson_residuals.pdf")   # open device
acf(df_model$pearson_res, 
    main="ACF of Pearson Residuals", 
    lag.max=30)                    # adjust lag.max as needed
dev.off()                          # close device

# 7.3 Observed vs. Predicted Counts

# Expected count under the hurdle-Poisson model
df_model$expected_hat <-
  (1 - df_model$pi_hat) *
  df_model$mu_hat /
  (-expm1(-df_model$mu_hat))

# Plot observed vs. predicted
ggplot(df_model, aes(x=date)) +
  geom_line(aes(y=deaths, color="Observed")) +
  geom_line(aes(y=expected_hat, color="Predicted")) +
  labs(title="Observed vs Predicted COVID-19 Deaths",
       y="Deaths", x="Date", color="Legend") +
  theme_minimal()

# Save observed vs predicted plot
obs_pred_plot <- ggplot(df_model, aes(x=date)) +
  geom_line(aes(y=deaths, color="Observed")) +
  geom_line(aes(y=expected_hat, color="Predicted")) +
  labs(title="Observed vs Predicted COVID-19 Deaths",
       y="Deaths", x="Date", color="Legend") +
  theme_minimal()

ggsave("observed_pred_plot.pdf", plot = obs_pred_plot, width=7, height=5)


#7.4 Zero-inflation Check

obs_zero <- sum(df_model$deaths == 0)
pred_zero <- sum(df_model$pi_hat)

cat("Observed zeros:", obs_zero, "\n")
cat("Predicted zeros (sum of pi_t):", round(pred_zero,2), "\n")

# Create a bar plot to compare observed vs predicted zeros

zero_plot <- ggplot(data.frame(
  Type = c("Observed Zeros", "Predicted Zeros"),
  Count = c(obs_zero, round(pred_zero,2))
), aes(x=Type, y=Count, fill=Type)) +
  geom_col(width=0.6) +
  coord_cartesian(ylim = c(0, 300)) +   # ✅ UPDATED SCALE
  geom_text(aes(label=round(Count,1)), vjust=-0.5) +
  labs(title="Observed vs Predicted Zeros",
       y="Count of Zeros", x="") +
  theme_minimal() +
  theme( legend.position = "none",
         panel.grid = element_blank())

ggsave("zero_fit_plot.pdf", plot = zero_plot, width=6, height=4)

nrow(df_model)

print(xtable(data.frame(
  Observed_Zeros = obs_zero,
  Predicted_Zeros = round(pred_zero, 2)
)), include.rownames = FALSE)

# 7.5 PIT or Uniformity Check (Optional)

# Randomized quantile residuals
library(pscl)

# Create function for Hurdle-Poisson residuals
rand_quant_resid <- function(deaths, mu_hat, pi_hat){
  n <- length(deaths)
  r <- numeric(n)
  for(i in 1:n){
    u <- runif(1)
    if(deaths[i] == 0){
      r[i] <- qnorm(u * pi_hat[i])
    } else {
      Fyt <- ppois(deaths[i]-1, lambda=mu_hat[i])
      r[i] <- qnorm(pi_hat[i] + (1 - pi_hat[i]) * Fyt + u*(1 - pi_hat[i])*dpois(deaths[i], mu_hat[i]))
    }
  }
  return(r)
}

df_model$rand_resid <- rand_quant_resid(df_model$deaths, df_model$mu_hat, df_model$pi_hat)

# --- 7.6 Randomized quantile residual diagnostics (lean set) ---

# 1) Histogram
library(ggplot2)
hist_res <- ggplot(df_model, aes(x = rand_resid)) +
  geom_histogram(binwidth = 0.5, color = "black", fill = "grey70") +
  labs(title = "Randomized Quantile Residuals",
       x = "Residual", y = "Frequency") +
  theme_minimal(base_size = 10)
ggsave("rand_resid_hist.pdf", plot = hist_res, width = 6.5, height = 4.5)

# 2) QQ plot
qq_res <- ggplot(df_model, aes(sample = rand_resid)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Randomized Quantile Residuals",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal(base_size = 10)
ggsave("rand_resid_qq.pdf", plot = qq_res, width = 6.5, height = 4.5)


################################################################################
# Monte Carlo: MLE validation for AR(1) Hurdle-Poisson
# - Uses your fitted 'est' as baseline true parameters
# - Sweeps sample sizes (N = 20,50,100,200,500)
# - For each sample size uses 3 true-parameter vectors:
#     1) est (fitted)
#     2) est - 0.05
#     3) est + 0.05
# - Records: Bias, RMSE, MAD, Coverage (95% CI from numeric Hessian), success rate
# - Exports CSV and LaTeX tables per sample size
################################################################################

library(dplyr)
library(numDeriv)
library(xtable)


# ---------- User configurable ----------
n_sim <- 500               # number of Monte Carlo replications (adjust if long)
sample_sizes <- c(20, 50, 100, 200, 500)
perturb <- 0.05            # amount to subtract/add from fitted est to create 3 param sets
out_dir <- "."             # save files here
set.seed(12345)
# ---------------------------------------

# Check that est and start_par and neg_logLik_AR1_HP exist
if(!exists("est") || !exists("start_par") || !exists("neg_logLik_AR1_HP")){
  stop("Please run the model-fitting code first so 'est', 'start_par', and 'neg_logLik_AR1_HP' exist.")
}


# Ensure est is named and contains required parameter names: alpha,eta,beta,gamma,theta,psi
req_names <- c("alpha","eta","beta","gamma","theta","psi")
if(!all(req_names %in% names(est))){
  stop("The fitted 'est' must contain named entries: alpha, eta, beta, gamma, theta, psi")
}


# Create three true-parameter vectors: est, est - perturb, est + perturb
base_params <- est[req_names]
param_sets <- list(
  "base" = base_params,
  "minus" = base_params - perturb,
  "plus"  = base_params + perturb
)

# Keep names
param_sets <- lapply(param_sets, function(x){ names(x) <- req_names; x })

# Helper: simulate one time-series from AR(1) HP given parameter vector (length n_time)
simulate_ar1_hp <- function(par, n_time){
  # par: named numeric vector alpha,eta,beta,gamma,theta,psi (psi on real scale)
  alpha <- as.numeric(par["alpha"]); eta <- as.numeric(par["eta"])
  beta  <- as.numeric(par["beta"]);  gamma<- as.numeric(par["gamma"])
  theta <- as.numeric(par["theta"]); psi  <- as.numeric(par["psi"])
  
  sim <- numeric(n_time)
  pi_prev <-  mean(df_model$deaths == 0) # initialize from data empirical
  mu_prev <-  mean(df_model$deaths[df_model$deaths>0])
  if(is.na(mu_prev) || mu_prev <= 0) mu_prev <- 1.0
  
  for(t in 1:n_time){
    x_lag <- if(t==1) 0 else sim[t-1]
    # zero probability
    linear_pi <- alpha + eta * pi_prev + beta * x_lag
    linear_pi <- pmin(pmax(linear_pi, -30), 30)
    pi_t <- 1/(1+exp(-linear_pi))
    
    # positive mean
    linear_mu <- gamma + theta * log(mu_prev + 1e-8) + psi * x_lag
    linear_mu <- pmin(pmax(linear_mu, -30), 30)
    mu_t <- exp(linear_mu)
    
    # hurdle draw
    if(runif(1) < pi_t){
      sim[t] <- 0
    } else {
      # Zero-truncated Poisson: positive component must be > 0
      repeat {
        y_pos <- rpois(1, lambda = mu_t)
        if(y_pos > 0) break
      }
      sim[t] <- y_pos
    }
    pi_prev <- pi_t; mu_prev <- mu_t
  }
  return(sim)
}


# Storage and results loop
for(N in sample_sizes){
  message("Running simulations for sample size N = ", N)
  
  
  # --- Checkpoint: skip if combined file already exists ---
  checkpoint_file <- file.path(out_dir, paste0("simulation_results_combined_N", N, ".tex"))
  if (file.exists(checkpoint_file)) {
    message("  Skipping N = ", N, " because results already exist: ", checkpoint_file)
    next
  }
  
  
  
  # Prepare containers: store estimates and SEs per sim and per param_set
  all_results <- list()
  
  for(ps_name in names(param_sets)){
    true_par <- param_sets[[ps_name]]
    message(" Parameter set: ", ps_name, " (", paste(round(true_par,4), collapse=", "), ")")
    
    # matrix to collect estimates, se, converged
    est_mat <- matrix(NA, nrow = n_sim, ncol = length(req_names))
    se_mat  <- matrix(NA, nrow = n_sim, ncol = length(req_names))
    conv_vec <- logical(n_sim)
    
    colnames(est_mat) <- req_names
    colnames(se_mat)  <- req_names
    
    pb <- txtProgressBar(min=0, max=n_sim, style=3)
    for(sim in 1:n_sim){
      # 1) simulate data of length N
      sim_series <- simulate_ar1_hp(true_par, n_time = N)
      sim_df <- data.frame(deaths = sim_series,
                           deaths_lag1 = dplyr::lag(sim_series, 1)) %>% na.omit()
      
      # 2) fit using nlminb starting at start_par (we reuse same start_par structure)
      fit_try <- tryCatch({
        out <- nlminb(start = start_par,
                      objective = function(p) neg_logLik_AR1_HP(p, sim_df),
                      control = list(iter.max = 20000, eval.max = 20000))
        out
      }, error = function(e) {
        return(NULL)
      })
      
      if(is.null(fit_try) || !is.numeric(fit_try$objective)){
        # failed fit => record NA
        conv_vec[sim] <- FALSE
      } else {
        # extract par (it might be named or not)
        par_hat <- if(!is.null(fit_try$par)) fit_try$par else if(!is.null(fit_try$estimate)) fit_try$estimate else NULL
        if(is.null(par_hat)){
          conv_vec[sim] <- FALSE
        } else {
          # Ensure names and ordering
          if(is.null(names(par_hat))) names(par_hat) <- names(start_par)
          
          # Keep only required params (alpha,eta,beta,gamma,theta,psi)
          par_hat <- par_hat[req_names]
          
          # numeric Hessian to estimate vcov (try / catch)
          hess_try <- tryCatch({
            numDeriv::hessian(func = function(p) neg_logLik_AR1_HP(p, sim_df), x = par_hat)
          }, error = function(e) NULL)
          
          if(!is.null(hess_try) && all(is.finite(hess_try))){
            vcov_try <- tryCatch({ solve(hess_try) }, error = function(e) NULL)
            if(!is.null(vcov_try) && all(is.finite(diag(vcov_try)))){
              se_hat <- sqrt(pmax(0, diag(vcov_try)))
            } else {
              se_hat <- rep(NA, length(par_hat))
            }
          } else {
            se_hat <- rep(NA, length(par_hat))
          }
          
          est_mat[sim, ] <- as.numeric(par_hat)
          se_mat[sim, ]  <- se_hat
          conv_vec[sim]  <- TRUE
        }
      }
      setTxtProgressBar(pb, sim)
    } # end sims
    close(pb)
    
    # Post-process: drop failed fits
    ok_idx <- which(conv_vec & apply(est_mat,1,function(x) all(!is.na(x))))
    n_ok <- length(ok_idx)
    success_rate <- n_ok / n_sim
    
    if(n_ok == 0){
      warning("No successful fits for N=", N, " param set ", ps_name)
      # create NA result table
      summary_df <- data.frame(Parameter = req_names,
                               True = as.numeric(true_par),
                               Bias = NA, RMSE = NA, MAD = NA, Coverage = NA)
    } else {
      est_ok <- est_mat[ok_idx, , drop=FALSE]
      se_ok  <- se_mat[ok_idx, , drop=FALSE]
      
      # convert to matrices (n_ok x p)
      # Compute Bias, RMSE, MAD (across successful fits)
      mean_est <- colMeans(est_ok, na.rm = TRUE)
      bias <- mean_est - as.numeric(true_par)
      rmse <- sqrt(colMeans((t(t(est_ok) - as.numeric(true_par)))^2, na.rm = TRUE))
      mad  <- colMeans(abs(t(t(est_ok) - as.numeric(true_par))), na.rm = TRUE)
      
      # Coverage: for each sim compute 95% CI using se from that sim (if se available)
      coverage_vec <- numeric(length(req_names))
      for(j in seq_along(req_names)){
        ses_j <- se_ok[, j]
        ests_j <- est_ok[, j]
        valid_j <- which(!is.na(ses_j) & ses_j > 0 & !is.na(ests_j))
        if(length(valid_j) > 0){
          lower <- ests_j[valid_j] - 1.96 * ses_j[valid_j]
          upper <- ests_j[valid_j] + 1.96 * ses_j[valid_j]
          coverage_vec[j] <- mean( (as.numeric(true_par[j]) >= lower) & (as.numeric(true_par[j]) <= upper) )
        } else {
          coverage_vec[j] <- NA
        }
      }
      
      summary_df <- data.frame(Parameter = req_names,
                               True = round(as.numeric(true_par), 5),
                               Bias = round(bias, 4),
                               RMSE = round(rmse, 4),
                               MAD = round(mad, 4),
                               Coverage = round(coverage_vec, 3),
                               SuccessRate = round(success_rate, 3),
                               stringsAsFactors = FALSE)
    }
    
    # Save per param-set table and whole raw estimates
    fname_csv <- file.path(out_dir, paste0("sim_raw_N", N, "_", ps_name, ".csv"))
    write.csv(data.frame(est_mat), file = fname_csv, row.names = FALSE)
    
    fname_summary_csv <- file.path(out_dir, paste0("sim_summary_N", N, "_", ps_name, ".csv"))
    write.csv(summary_df, file = fname_summary_csv, row.names = FALSE)
    
    # Save LaTeX table for this (N, param set)
    xt_sim <- xtable(summary_df,
                     caption = paste0("Simulation results (N=", N, ", set=", ps_name, "): Bias, RMSE, MAD, Coverage"),
                     label = paste0("tab:sim_N", N, "_", ps_name))
    print(xt_sim, include.rownames = FALSE,
          file = file.path(out_dir, paste0("simulation_results_N", N, "_", ps_name, ".tex")))
    
    all_results[[ps_name]] <- list(summary = summary_df,
                                   rawfile = fname_csv,
                                   summaryfile = fname_summary_csv)
  } # end param_sets
  
  # Combine the three parameter-set summaries into one LaTeX file for convenience
  combined <- do.call(rbind, lapply(names(all_results), function(nm){
    data.frame(Set = nm, all_results[[nm]]$summary, stringsAsFactors = FALSE)
  }))
  write.csv(combined, file = file.path(out_dir, paste0("sim_combined_N", N, ".csv")), row.names = FALSE)
  print(xtable(combined,
               caption = paste0("Combined sim results for N=", N),
               label = paste0("tab:sim_combined_N", N)),
        include.rownames = FALSE,
        file = file.path(out_dir, paste0("simulation_results_combined_N", N, ".tex")))
  
} # end sample_sizes

message("Simulation finished. Files written to: ", normalizePath(out_dir))

################################################################################
# Step 1: Load packages
################################################################################
library(dplyr)
library(ggplot2)
library(readr)

################################################################################
# Step 2: Read all combined simulation results
################################################################################
# Adjust pattern if your filenames differ
files <- list.files(pattern = "sim_combined_N[0-9]+\\.csv$")

# Extract N from filename
extract_N <- function(fname){
  as.numeric(gsub("[^0-9]", "", fname))
}

all_sims <- do.call(rbind, lapply(files, function(f){
  N_val <- extract_N(f)
  df <- read_csv(f, show_col_types = FALSE)
  df$N <- N_val
  return(df)
}))

# Ensure ordering
all_sims <- all_sims %>% arrange(N, Set, Parameter)

################################################################################
# Step 3: Quick sanity check
################################################################################
print(head(all_sims))
table(all_sims$N, all_sims$Set)

################################################################################
# Step 4: Make plots for Bias, RMSE, Coverage, SuccessRate
################################################################################
# Bias vs. N
p_bias <- ggplot(all_sims, aes(x = N, y = Bias, color = Parameter, linetype = Set)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "Bias vs Sample Size",
       y = "Bias")

print(p_bias)  # shows on screen
ggsave("plot_bias.png", plot = p_bias, width = 7, height = 5, dpi = 300)  # SAVE PLOT

# RMSE vs. N
p_rmse <- ggplot(all_sims, aes(x = N, y = RMSE, color = Parameter, linetype = Set)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "RMSE vs Sample Size",
       y = "RMSE")

print(p_rmse)
ggsave("plot_rmse.png", plot = p_rmse, width = 7, height = 5, dpi = 300)

# Coverage vs. N
p_cov <- ggplot(all_sims, aes(x = N, y = Coverage, color = Parameter, linetype = Set)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Coverage vs Sample Size",
       y = "Coverage Probability")

print(p_cov)
ggsave("plot_coverage.png", plot = p_cov, width = 7, height = 5, dpi = 300)

# Success Rate vs. N
p_succ <- ggplot(all_sims %>% distinct(N, Set, SuccessRate, .keep_all = TRUE),
                 aes(x = N, y = SuccessRate, color = Set)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "Convergence Success Rate vs Sample Size",
       y = "Success Rate")

print(p_succ)
ggsave("plot_successrate.png", plot = p_succ, width = 7, height = 5, dpi = 300)


# --- Packages ---
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(xtable)

# --- Function to parse a LaTeX table into a clean data frame ---
read_sim_tex <- function(file, N) {
  raw <- read_lines(file)
  
  # Keep only rows with actual data (contain "&" and not \hline, not caption)
  data_rows <- raw[grepl("&", raw) & !grepl("hline", raw) & !grepl("caption", raw)]
  
  # Split by "&" and trim spaces
  df <- str_split(data_rows, "&") %>% 
    lapply(str_trim) %>% 
    do.call(rbind, .) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  
  # First two cols are "Set" and "Parameter"
  # Adjust depending on your tex structure (from your pasted N=500 table)
  colnames(df) <- c("Set", "Parameter", "True", "Bias", "RMSE", "MAD", "Coverage", "SuccessRate")
  
  # Convert numeric cols
  df <- df %>%
    mutate(across(True:Coverage, as.numeric)) %>%      # convert other numeric cols normally
    mutate(SuccessRate = as.numeric(gsub("[^0-9.-]", "", SuccessRate))) %>% # clean SuccessRate
    mutate(N = N)
  
  return(df)
}

# --- Read all N files ---
Ns <- c(20, 50, 100, 200, 500)
files <- paste0("simulation_results_combined_N", Ns, ".tex")
all_results <- bind_rows(mapply(read_sim_tex, files, Ns, SIMPLIFY = FALSE))

export_sim_table <- function(df, set_name, file_name, caption_label) {
  
  # Compute Estimate column
  df <- df %>% mutate(Estimate = True + Bias)
  
  ################################################################################
  # Step 5: Export LaTeX Table of Simulation Results
  ################################################################################
  library(xtable)
  library(tidyr)
  ################################################################################
  
  # Summarise across parameters and sample sizes
  
  df$Parameter <- recode(df$Parameter,
                         alpha = "$\\alpha$",
                         eta   = "$\\eta$",
                         beta  = "$\\beta$",
                         gamma = "$\\gamma$",
                         theta = "$\\theta$",
                         psi   = "$\\psi$")
  
  
  # Pivot to wide format
  latex_tbl <- df %>%
    select(-Set, -MAD) %>%
    pivot_longer(cols = c(True, Estimate, Bias, RMSE, Coverage, SuccessRate),
                 names_to = "Metric", values_to = "Value") %>%
    pivot_wider(names_from = N, values_from = Value) %>%
    arrange(Parameter, factor(Metric,
                              levels = c("True", "Estimate", "Bias", "RMSE", "Coverage", "SuccessRate"))) %>%
    mutate(Metric = recode(Metric,
                           True = "\\textbf{True}",
                           Estimate = "\\textbf{Estimate}",
                           Bias = "\\textbf{Bias}",
                           RMSE = "\\textbf{RMSE}",
                           Coverage = "\\textbf{Coverage}",
                           SuccessRate = "\\textbf{Success Rate}"))
  
  # Add N labels to columns
  colnames(latex_tbl)[-(1:2)] <- paste0("N(", colnames(latex_tbl)[-(1:2)], ")")
  
  # Blank repeated parameter names
  latex_tbl$Parameter <- ave(latex_tbl$Parameter, latex_tbl$Parameter,
                             FUN = function(x) c(x[1], rep("", length(x)-1)))
  
  # Add midrules between parameters
  param_breaks <- which(latex_tbl$Parameter != "")[-1] - 1
  param_breaks <- c(param_breaks, nrow(latex_tbl))
  addtorow <- list(
    pos = as.list(param_breaks),
    command = rep("\\midrule\n", length(param_breaks))
  )
  
  # Export LaTeX table
  print(
    xtable(latex_tbl,
           caption = caption_label,
           label = paste0("tab:sim_results_", set_name),
           digits = 3),
    include.rownames = FALSE,
    floating = FALSE,
    booktabs = TRUE,
    sanitize.text.function = identity,
    add.to.row = addtorow,
    file = file_name
  )
}

# Base set
base_results <- all_results %>% filter(Set == "base")
export_sim_table(base_results, "base", "sim_results_base.tex",
                 "Simulation results for the AR(1) Hurdle Poisson Model. True values from the base set.")

# Minus set
minus_results <- all_results %>% filter(Set == "minus")
export_sim_table(minus_results, "minus", "sim_results_minus.tex",
                 "Simulation results for the AR(1) Hurdle Poisson Model. True values from the minus set.")

# Plus set
plus_results <- all_results %>% filter(Set == "plus")
export_sim_table(plus_results, "plus", "sim_results_plus.tex",
                 "Simulation results for the AR(1) Hurdle Poisson Model. True values from the plus set.")



# --- Split into training / hold-out ---
n <- nrow(df_model)
train_df <- df_model[1:(n-14), ]
test_df  <- df_model[(n-13):n, ]

# --- Use confirmed model estimates for forecasting ---
fitted_params <- est

print(fitted_params)

# --- Forecast function ---
forecast_HP_AR1_mc <- function(fit_params, last_observed, horizon, sims = 500) {
  all_sims <- matrix(NA, nrow = sims, ncol = horizon)
  
  for (s in 1:sims) {
    pi_prev <- ifelse(last_observed == 0, 1, 0)
    mu_prev <- max(last_observed, 1)
    last_obs <- last_observed
    
    for (t in 1:horizon) {
      # zero prob
      linear_pi <- fit_params["alpha"] + fit_params["eta"] * pi_prev + fit_params["beta"] * last_obs
      pi_t <- 1 / (1 + exp(-linear_pi))
      
      # positive mean
      linear_mu <- fit_params["gamma"] + fit_params["theta"] * log(mu_prev + 1e-8) + fit_params["psi"] * last_obs
      mu_t <- exp(linear_mu)
      
      # simulate from the hurdle distribution
      y_t <- if(runif(1) < pi_t){
        y_t <- 0
      } else {
        repeat {
          y_pos <- rpois(1, mu_t)
          if(y_pos > 0) break
        }
        y_t <- y_pos
      }
      
      all_sims[s, t] <- y_t
      
      # update states
      pi_prev <- pi_t
      mu_prev <- mu_t
      last_obs <- y_t
    }
  }
  
  mean_forecast <- apply(all_sims, 2, mean)
  lower <- apply(all_sims, 2, quantile, probs = 0.025)
  upper <- apply(all_sims, 2, quantile, probs = 0.975)
  
  return(data.frame(
    DayAhead = 1:horizon,
    Mean = round(mean_forecast, 1),
    Lower = round(lower, 1),
    Upper = round(upper, 1)
  ))
}

# --- Generate forecasts ---
last_obs <- tail(train_df$deaths, 1)
fc_7  <- forecast_HP_AR1_mc(fitted_params, last_obs, horizon = 7, sims = 500)
fc_14 <- forecast_HP_AR1_mc(fitted_params, last_obs, horizon = 14, sims = 500)

# Attach actuals
fc_7$Actual  <- test_df$deaths[1:7]
fc_14$Actual <- test_df$deaths[1:14]


# --- Quick check-----
print(fc_7)
print(fc_14)


# --- Export to LaTeX ---

# Save 7-day forecast

capture.output(
  print(
    xtable(fc_7),
    include.rownames = FALSE,
    booktabs = TRUE,
    floating = FALSE,         # <---- IMPORTANT, prevents \begin{table}
    tabular.environment = "tabular",
    hline.after = c(-1, 0, nrow(fc_7)) # top, header, bottom lines
  ),
  file = "forecast_7day.tex"
)

# Save 14-day forecast
capture.output(
  print(
    xtable(fc_14),
    include.rownames = FALSE,
    booktabs = TRUE,
    floating = FALSE,         # <---- prevents extra table wrapper
    tabular.environment = "tabular",
    hline.after = c(-1, 0, nrow(fc_14))
  ),
  file = "forecast_14day.tex"
)



# --- Load packages ---
library(dplyr)
library(readr)
library(stringr)
library(xtable)

# --- Function to clean and extract numbers from .tex forecast files ---
read_forecast_tex <- function(file) {
  lines <- readLines(file)
  
  # Keep only rows with '&'
  rows <- lines[grepl("&", lines)]
  
  # Split and clean
  data <- str_split(rows, "&") %>%
    lapply(function(x) str_trim(gsub("\\\\.*", "", x))) %>%
    do.call(rbind, .)
  
  df <- as.data.frame(data, stringsAsFactors = FALSE)
  
  # Convert numeric columns
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
  
  # Correct column names: DayAhead | Forecast | Lower | Upper | Observed
  colnames(df) <- c("Day", "Forecast", "Lower", "Upper", "Observed")
  
  return(df)
}

# --- Function to calculate forecast accuracy metrics ---
forecast_metrics <- function(df) {
  obs   <- df$Observed
  pred  <- df$Forecast
  lower <- df$Lower
  upper <- df$Upper
  
  valid_idx <- !is.na(obs) & !is.na(pred)
  obs   <- obs[valid_idx]
  pred  <- pred[valid_idx]
  lower <- lower[valid_idx]
  upper <- upper[valid_idx]
  
  # RMSE & MAE
  RMSE <- sqrt(mean((pred - obs)^2, na.rm = TRUE))
  MAE  <- mean(abs(pred - obs), na.rm = TRUE)
  
  # MAPE (skip zeros)
  nonzero <- obs != 0
  if (any(nonzero)) {
    MAPE <- mean(abs((pred[nonzero] - obs[nonzero]) / obs[nonzero]), na.rm = TRUE) * 100
  } else {
    MAPE <- NA
  }
  
  # sMAPE (handle 0/0 safely)
  SMAPE <- mean(
    ifelse((abs(obs) + abs(pred)) == 0, 0,
           abs(pred - obs) / ((abs(obs) + abs(pred)) / 2)),
    na.rm = TRUE
  ) * 100
  
  # Coverage
  Coverage <- mean(obs >= lower & obs <= upper, na.rm = TRUE) * 100
  
  return(data.frame(
    RMSE = RMSE,
    MAE = MAE,
    MAPE = MAPE,
    SMAPE = SMAPE,
    Coverage = Coverage
  ))
}

# --- Read forecasts ---
fcst7  <- read_forecast_tex("forecast_7day.tex")
fcst14 <- read_forecast_tex("forecast_14day.tex")

# --- Compute metrics ---
metrics7  <- forecast_metrics(fcst7)
metrics14 <- forecast_metrics(fcst14)

metrics_all <- rbind(
  cbind(Horizon = "7-day", metrics7),
  cbind(Horizon = "14-day", metrics14)
)

# --- Export LaTeX table safely ---
# Avoids "extra table environment" by using only tabular environment
print(
  xtable(metrics_all,
         digits = c(0,0,3,3,2,2,2)),
  include.rownames = FALSE,
  tabular.environment = "tabular",
  hline.after = c(-1,0,nrow(metrics_all)),  # adds top, mid, bottom rules
  file = "forecast_metrics.tex"
)


print(metrics_all)

####################CREATING FORECAST PLOT


# --- 7-DAY FORECAST PLOT ---

pdf("forecast_7day_plot.pdf", width = 9, height = 5)

# Increase right margin and allow drawing outside plot
par(mar = c(5, 5, 4, 10), xpd = TRUE)

fc_7$Time <- 1:nrow(fc_7)

# Plot observed
plot(fc_7$Time, fc_7$Actual,
     type = "l", lwd = 2, col = "black",
     xlab = "Days Ahead",
     ylab = "Deaths",
     ylim = range(c(fc_7$Actual, fc_7$Lower, fc_7$Upper)),
     main = "7-Day Forecast vs Observed")

# Add prediction interval shading (very important)
polygon(c(fc_7$Time, rev(fc_7$Time)),
        c(fc_7$Lower, rev(fc_7$Upper)),
        col = rgb(1, 0, 0, 0.2), border = NA)

# Add forecast mean
lines(fc_7$Time, fc_7$Mean, lty = 2, lwd = 2, col = "blue")

# Add interval boundaries (bold)
lines(fc_7$Time, fc_7$Lower, col = "red", lwd = 3)
lines(fc_7$Time, fc_7$Upper, col = "red", lwd = 3)

# Add observed again (to keep it on top)
lines(fc_7$Time, fc_7$Actual, lwd = 2, col = "black")

# Legend outside plot
legend("topright",
       inset = c(-0.35, 0),
       legend = c("Observed", "Forecast Mean", "Prediction Interval"),
       col = c("black", "blue", "red"),
       lty = c(1, 2, 1),
       lwd = c(2, 2, 3),
       bty = "n")

dev.off()

# --- 14-DAY FORECAST PLOT ---


pdf("forecast_14day_plot.pdf", width = 9, height = 5)

par(mar = c(5, 5, 4, 10), xpd = TRUE)

fc_14$Time <- 1:nrow(fc_14)

plot(fc_14$Time, fc_14$Actual,
     type = "l", lwd = 2, col = "black",
     xlab = "Days Ahead",
     ylab = "Deaths",
     ylim = range(c(fc_14$Actual, fc_14$Lower, fc_14$Upper)),
     main = "14-Day Forecast vs Observed")

# Shaded interval
polygon(c(fc_14$Time, rev(fc_14$Time)),
        c(fc_14$Lower, rev(fc_14$Upper)),
        col = rgb(1, 0, 0, 0.2), border = NA)

# Forecast mean
lines(fc_14$Time, fc_14$Mean, lty = 2, lwd = 2, col = "blue")

# Interval bounds
lines(fc_14$Time, fc_14$Lower, col = "red", lwd = 3)
lines(fc_14$Time, fc_14$Upper, col = "red", lwd = 3)

# Observed (on top)
lines(fc_14$Time, fc_14$Actual, lwd = 2, col = "black")

# Legend outside
legend("topright",
       inset = c(-0.35, 0),
       legend = c("Observed", "Forecast Mean", "Prediction Interval"),
       col = c("black", "blue", "red"),
       lty = c(1, 2, 1),
       lwd = c(2, 2, 3),
       bty = "n")

dev.off()


