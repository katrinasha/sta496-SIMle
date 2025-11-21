# High-dimensional Gaussian approximation demo (R)
# Simulates locally dependent d-dimensional vectors X_i,
# compares distribution of max-component of normalized sum to Gaussian approx.

set.seed(2025)
library(MASS)   # for mvrnorm
library(ggplot2)

# Parameters
n <- 200        # number of summands (observations)
d <- 100        # dimension of each vector (high-d)
rho <- 0.4      # local dependence strength across i (one-step)
nrep <- 2000    # number of Monte Carlo experiments to estimate empirical distribution
ngsim <- 2000   # number of Gaussian draws for Gaussian approx

# Function to generate one sequence of X_i (i=1..n) with local dependence across i
# Each epsilon_i is i.i.d. N(0, I_d); X_i = epsilon_i + rho * epsilon_{i-1}
generate_X_sequence <- function(n, d, rho) {
  eps <- matrix(rnorm((n+1)*d), nrow = n+1, ncol = d) # eps_0..eps_n
  X <- matrix(0, nrow = n, ncol = d)
  for (i in 1:n) {
    X[i, ] <- eps[i+1, ] + rho * eps[i, ]  # uses eps_{i} and eps_{i-1}
  }
  X
}

# Monte Carlo: compute empirical distribution of T = max_j W_j
Tvals <- numeric(nrep)
for (rep in 1:nrep) {
  X <- generate_X_sequence(n, d, rho)
  W <- colSums(X) / sqrt(n)        # normalized sum, length d
  Tvals[rep] <- max(W)
  if (rep %% 500 == 0) cat("done rep", rep, "\n")
}

# Estimate covariance Sigma of the sum W from one long replicate or theoretical:
# We'll estimate Sigma empirically by averaging many independent sequences (cheap approx).
# Alternatively, approximate Sigma by Var(colSums(X)/sqrt(n)) across additional sims.
mSigma_rep <- 1000
Slist <- matrix(0, nrow = mSigma_rep, ncol = d)
for (r in 1:mSigma_rep) {
  X <- generate_X_sequence(n, d, rho)
  Slist[r, ] <- colSums(X) / sqrt(n)
}
Sigma_hat <- cov(Slist)  # d x d estimated covariance of W

# Gaussian approximation: draw from N(0, Sigma_hat) and compute max coordinate
gmax <- numeric(ngsim)
cholS <- chol(Sigma_hat + 1e-8 * diag(d))  # add tiny ridge for numerical stability
for (k in 1:ngsim) {
  z <- rnorm(d)
  g <- cholS %*% z
  gmax[k] <- max(as.numeric(g))
}

# Compare empirical CDFs for max statistic
ecdf_true <- ecdf(Tvals)
ecdf_gauss <- ecdf(gmax)

# Compute Kolmogorov distance between the two empirical CDFs at combined grid
grid <- sort(unique(c(Tvals, gmax)))
d_ks <- max(abs(ecdf_true(grid) - ecdf_gauss(grid)))
cat("Kolmogorov distance (empirical) between max(T) and Gaussian max:", round(d_ks,4), "\n")

# Plot the two CDFs
df_plot <- data.frame(
  value = c(grid, grid),
  cdf = c(ecdf_true(grid), ecdf_gauss(grid)),
  type = rep(c("True max (dependent sum)", "Gaussian approx"), each = length(grid))
)

ggplot(df_plot, aes(x = value, y = cdf, color = type)) +
  geom_line(size = 1) +
  labs(title = "CDF Comparison: max-coordinate of normalized sum vs Gaussian approx",
       subtitle = paste0("n=", n, ", d=", d, ", rho=", rho, ", nrep=", nrep),
       x = "value", y = "CDF") +
  theme_minimal(base_size = 13)

# Also show QQ-plot of quantiles
q_emp <- quantile(Tvals, probs = seq(0.01,0.99,by=0.01))
q_gau <- quantile(gmax, probs = seq(0.01,0.99,by=0.01))
qqdf <- data.frame(emp = q_emp, gau = q_gau)
ggplot(qqdf, aes(x = gau, y = emp)) +
  geom_point() + geom_abline(slope=1, intercept=0, linetype=2) +
  labs(title = "QQ-plot: Empirical max vs Gaussian max",
       subtitle = paste0("Kolmogorov dist = ", round(d_ks,4)),
       x = "Gaussian max quantiles", y = "Empirical max quantiles") +
  theme_minimal(base_size = 13)

# Print a short summary
cat("\nSummary statistics:\n")
cat("Empirical max: mean=", round(mean(Tvals),3), " sd=", round(sd(Tvals),3), "\n")
cat("Gaussian max:  mean=", round(mean(gmax),3), " sd=", round(sd(gmax),3), "\n")

##################### Example 2: L2 norm ###########################
set.seed(2025)
library(MASS)
library(ggplot2)

# Parameters
n <- 200   # sample size
d <- 50    # dimension
rho <- 0.4 # local dependence
nrep <- 2000
ngsim <- 2000

# Generate sequence with local dependence
generate_X_sequence <- function(n, d, rho) {
  eps <- matrix(rnorm((n+1)*d), nrow = n+1, ncol = d)
  X <- matrix(0, nrow = n, ncol = d)
  for (i in 1:n) {
    X[i, ] <- eps[i+1, ] + rho * eps[i, ]
  }
  X
}

# Empirical distribution of ||W||_2
Tvals <- numeric(nrep)
for (rep in 1:nrep) {
  X <- generate_X_sequence(n, d, rho)
  W <- colSums(X) / sqrt(n)
  Tvals[rep] <- sqrt(sum(W^2))  # L2 norm
}

# Estimate covariance Sigma
mSigma_rep <- 1000
Slist <- matrix(0, nrow = mSigma_rep, ncol = d)
for (r in 1:mSigma_rep) {
  X <- generate_X_sequence(n, d, rho)
  Slist[r, ] <- colSums(X) / sqrt(n)
}
Sigma_hat <- cov(Slist)

# Gaussian approximation: ||N(0, Sigma_hat)||_2
gvals <- numeric(ngsim)
cholS <- chol(Sigma_hat + 1e-8 * diag(d))
for (k in 1:ngsim) {
  z <- rnorm(d)
  g <- cholS %*% z
  gvals[k] <- sqrt(sum(g^2))
}

# Compare CDFs
ecdf_true <- ecdf(Tvals)
ecdf_gauss <- ecdf(gvals)
grid <- sort(unique(c(Tvals, gvals)))
d_ks <- max(abs(ecdf_true(grid) - ecdf_gauss(grid)))
cat("Kolmogorov distance between empirical L2 norm and Gaussian L2 norm:", round(d_ks,4), "\n")

# Plot
df_plot <- data.frame(
  value = c(grid, grid),
  cdf = c(ecdf_true(grid), ecdf_gauss(grid)),
  type = rep(c("True L2 norm", "Gaussian approx"), each = length(grid))
)

ggplot(df_plot, aes(x = value, y = cdf, color = type)) +
  geom_line(size = 1) +
  labs(title = "Gaussian Approximation for L2 Norm of High-dim Sum",
       subtitle = paste0("n=", n, ", d=", d, ", rho=", rho),
       x = "value", y = "CDF") +
  theme_minimal(base_size = 13)

################## Counterexample: median ################################
set.seed(2025)
library(MASS)

# Parameters
n <- 200
d <- 50
rho <- 0.4
nrep <- 2000
ngsim <- 2000

# Local dependence generator
generate_X_sequence <- function(n, d, rho) {
  eps <- matrix(rnorm((n+1)*d), nrow = n+1, ncol = d)
  X <- matrix(0, nrow = n, ncol = d)
  for (i in 1:n) {
    X[i, ] <- eps[i+1, ] + rho * eps[i, ]
  }
  X
}

# Empirical distribution of median of coordinates
Tvals <- numeric(nrep)
for (rep in 1:nrep) {
  X <- generate_X_sequence(n, d, rho)
  W <- colSums(X) / sqrt(n)
  Tvals[rep] <- median(W)  # median statistic
}

# Estimate covariance Sigma
mSigma_rep <- 1000
Slist <- matrix(0, nrow = mSigma_rep, ncol = d)
for (r in 1:mSigma_rep) {
  X <- generate_X_sequence(n, d, rho)
  Slist[r, ] <- colSums(X) / sqrt(n)
}
Sigma_hat <- cov(Slist)

# Gaussian approximation: median of N(0,Sigma_hat)
gvals <- numeric(ngsim)
cholS <- chol(Sigma_hat + 1e-8 * diag(d))
for (k in 1:ngsim) {
  z <- rnorm(d)
  g <- cholS %*% z
  gvals[k] <- median(as.numeric(g))
}

# Compare means and sds
cat("Empirical median: mean =", round(mean(Tvals),3), " sd =", round(sd(Tvals),3), "\n")
cat("Gaussian median:  mean =", round(mean(gvals),3), " sd =", round(sd(gvals),3), "\n")

# KS test
ks <- suppressWarnings(ks.test(Tvals, gvals))
cat("Kolmogorov-Smirnov distance:", ks$statistic, "\n")

# Plot densities
hist(Tvals, probability=TRUE, col=rgb(1,0,0,0.3), breaks=40,
     main="Median of Coordinates: True vs Gaussian approx",
     xlab="Median statistic")
hist(gvals, probability=TRUE, col=rgb(0,0,1,0.3), breaks=40, add=TRUE)
legend("topright", legend=c("Empirical", "Gaussian approx"),
       fill=c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)))

#################### Counterexample: non-local dependence of summands #########
set.seed(1)
library(MASS)

n    <- 300             # time points
p    <- 50              # dimension
nrep <- 3000            # Monte Carlo reps

b <- rep(1/sqrt(p), p)  # loading
a <- b / sqrt(sum(b^2)) # test direction (convex: linear functional)

Tvals <- numeric(nrep)

for (r in 1:nrep) {
  U  <- rt(1, df = 3)                 # common heavy-tailed factor (non-Gaussian)
  E  <- matrix(rnorm(n*p), n, p)      # i.i.d. Gaussian idiosyncratic noise
  X  <- sweep(E, 2, 0, "+") + matrix(U, n, 1) %*% t(b)  # X_i = b*U + ε_i
  W  <- colSums(X) / sqrt(n)          # normalized sum
  Tvals[r] <- sum(a * W)              # convex statistic: a^T W
}

# Gaussian proxy that matches mean and variance of T
mu  <- mean(Tvals); s2 <- var(Tvals)
Tgau <- rnorm(nrep, mean = 0, sd = sqrt(s2))

# KS distance; increase n and see it does NOT go to 0
suppressWarnings(ks.test(Tvals, Tgau))

