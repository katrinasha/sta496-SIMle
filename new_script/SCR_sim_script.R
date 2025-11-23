library(SIMle)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(parallel)
library(snow)

################### Error functions #####################
# tvAR2 helpers
a1_fun <- function(u) 0.3
a2_fun <- function(u) 0.3 * sin(2*pi*u)

# tvAR2
tvAR2_eps <- function(n, a1 = a1_fun, a2 = a2_fun, burn = 500L) {
  N <- n + burn
  eps <- rep(0, N)
  eta <- rnorm(N)
  for (t in 3:N) {
    u <- t/N
    eps[t] <- a1(u)*eps[t-1] + a2(u)*eps[t-2] + eta[t]
  }
  tail(eps, n) # remove the first 'burn' values to reduce the effect of initial values
}

# SETAR
setar_eps <- function(n, a1 = a1_fun, a2 = a2_fun, burn = 500L) {
  N <- n + burn
  eps <- rep(0, N)
  eta <- rnorm(N)
  for (t in 2:N) {
    u <- t/N
    a <- if (eps[t-1] >= 0) a1(u) else a2(u)
    eps[t] <- a*eps[t-1] + eta[t]
  }
  tail(eps, n)
}

# Bilinear
bilinear_eps <- function(n, a1 = a1_fun, a2 = a2_fun, burn = 500L) {
  N <- n + burn
  eps <- rep(0, N)
  eta <- rnorm(N)
  for (t in 2:N) {
    u <- t/N
    eps[t] <- (a1(u)*eta[t] + a2(u)) * eps[t-1] + eta[t]
  }
  tail(eps, n)
}

errors <- list(
  tvAR2 = tvAR2_eps,
  SETAR = setar_eps,
  bilinear = bilinear_eps
)

################### Main functions: m and sigma #####################
m1 <- function(t, x, delta) (1 + x^2)^(-4) + delta * sin(2*pi*t) * exp(-x^2)
s1 <- function(t, x)        1.5 * exp(-x^2/2) * (2 + sin(2*pi*t))
case1 <- list(m = m1, s = s1)

m2 <- function(t, x, delta) (delta * sin(2*pi*t) + 1) * exp(-x^2/2)
s2 <- function(t, x)        0.5 * exp(-x^2) * cos(2*pi*t) + 1
case2 <- list(m = m2, s = s2)


m3 <- function(t, x, delta) 2*t*(delta*exp(-2*t*x^2)) + pi^(-1/2)*exp(-x^2/2)
s3 <- function(t, x) {
  ax <- abs(x)
  if (ax <= 1)     return(0.7 * (1 + x^2))
  if (ax > 1 & 0 <= t & t < 0.5)   return(1.4) 
  return(2.0)
}
case3 <- list(m = m3, s = s3)

model_cases <- list(case1, case2, case3)


################ Simulate TS x_t #####################
sim_x <- function(n, delta=1, m_case, eps_case, x0=0, burn=500L) {
    stopifnot(m_case %in% 1:3)
    case_function <- model_cases[[m_case]]
    m_fun <- case_function$m
    s_fun <- case_function$s
    eps_fun <- errors[[eps_case]]

    # Generate x_t
    eps <- eps_fun(n, burn = burn)

    X <- numeric(n)
    x_prev <- x0
    data_used <- data.frame()
    for (t in 1:n) {
        u <- t/n
        m_val <- m_fun(u, x_prev, delta)
        s_var <- s_fun(u, x_prev)
        x_t <- m_val + s_var * eps[t]
        data_used <- rbind(data_used, data.frame(t=u, x = x_t, x_prev=x_prev, eps=eps[t]))
        X[t] <- x_t
        x_prev <- x_t

    }
    data_used %>% data.table()
    return(X)
}

# n=500
# m_case=1
# eps_case='tvAR2'
# set.seed(1)
# test_data <- sim_x(n, m_case = m_case, eps_case = eps_case)
# write.csv(test_data, file = 'method/test_data.csv', row.names = FALSE)
# plot(sim_x(n, m_case = m_case, eps_case = eps_case), type='l')

################### Verify SCR coverage #####################
extract_band <- function(df) {
  n_point <- length(df$x) / 3
  x_grid  <- df$x[1:n_point]
  idx_e   <- df$order == "esti"
  idx_u   <- df$order == "upper"
  idx_l   <- df$order == "lower"
  list(x = x_grid,
       esti  = df$y[idx_e],
       upper = df$y[idx_u],
       lower = df$y[idx_l])
}

auto_SCR_coverage <- function(n_rep, n, c, d, delta, b_time, b_timese, mp_type, 
                             type = 'fixt', ops = 'CV', per = 0.1, m_case, eps_case) {

    # Steps for checking coverage:
    #   (1) Define a grid for t on [0,1], ensure there is a point at each 0.01 interval
    #   (2) Move along t grid: For each t0 in the grid
    #       Estimate m(t0, x) and SCR bands
    #       Check if true m(t0, x) is within the bands for all x in the band grid
    #   (3) If yes for all t0, then mark this replication as covered
    #   (4) Repeat for n_rep replications, return the coverage proportion

    covered <- logical(n_rep)
    t_grid <- seq(0, 1, length.out = 101)
    for (s in 1:n_rep) {
        cat('Simulation ', s, ' / ', n_rep, '\n')

        set.seed(s)
        X <- sim_x(n, delta = delta, m_case = m_case, eps_case = eps_case, x0=0)
        ok_all_t <- TRUE

        for (t0 in t_grid) {
            cat('\tt0 = ', t0, '\n')

            res <- auto.SCR(X, c=c, d=d, b_time=b_time, b_timese=b_timese,
                            mp_type=mp_type, type='fixt', ops=ops, fix_num = t0, per=per)
            df <- res[[1]]
            band <- extract_band(df)
            m_true <- model_cases[[m_case]]$m
            m_vals <- m_true(t = rep(t0, length(band$x)), x = band$x, delta = delta)

            # Plot SCR and true m at this t0
            ggplot() +
                geom_line(aes(x = band$x, y = band$esti, color = 'Estimated')) +
                geom_line(aes(x = band$x, y = band$upper, color = 'Confidence region'), linetype = 'dashed') +
                geom_line(aes(x = band$x, y = band$lower, color = 'Confidence region'), linetype = 'dashed') +
                geom_line(aes(x = band$x, y = m_vals, color = 'True m')) +
                scale_color_manual(values = c('Estimated' = 'blue', 'Confidence region' = 'red', "True m" = 'green')) +
                ggtitle(paste0('Simulation ', s, ', t0=', round(t0, 3), 'c=', c, 'd=', d, 'mp_type=', mp_type)) +
                ylab('m(t0, x) and SCR') +
                theme(legend.position = 'top') + 
                labs(title = paste0('Simulation ', s, ', t0=', round(t0, 3), 'c=', c, 'd=', d, 'mp_type=', mp_type), 
                     subtitle = paste0('m_case=', m_case, ', eps_case=', eps_case))

            covered_t <- all(m_vals >= band$lower & m_vals <= band$upper)
            if (!covered_t) {
                ok_all_t <- FALSE
                break
            }
        }
        covered[s] <- ok_all_t
        cat('Current uncovered vs. covered = ', table(covered), '\n')
    }
    return(mean(covered))

}
##################### Example runs #####################
n_rep = 500
n = 500
c = 5 # consulted Professor Zhou for best choice for c and d
d = 5
delta = 1
b_time = 'Legen'
b_timese = 'Legen'
mp_type = 'algeb'
m_case = 1
eps_case = 'tvAR2'

# Issue: low empirical coverage
auto_SCR_coverage(n_rep = n_rep, n = n, c = c, d = d, delta = delta, b_time = b_time, b_timese = b_timese, 
                 mp_type = mp_type, m_case = m_case, eps_case = eps_case)

# Observed that SCR often fails to cover true m around x=0
#   See this issue for fixed t0=0 below
set.seed(2)
X <- sim_x(n, delta = delta, m_case = m_case, eps_case = eps_case, x0=0)
res_t <- auto.SCR(X, c=c, d=d, b_time=b_time, b_timese=b_timese,
                mp_type=mp_type, type='fixt', ops='CV', fix_num = 0, per=0.1)
df <- res_t[[1]]
band <- extract_band(df)
m_true <- model_cases[[m_case]]$m
m_vals <- m_true(t = rep(0, length(band$x)), x = band$x, delta = delta)
ggplot() +
    geom_line(aes(x = band$x, y = band$esti, color = 'Estimated')) +
    geom_line(aes(x = band$x, y = band$upper, color = 'Confidence region'), linetype = 'dashed') +
    geom_line(aes(x = band$x, y = band$lower, color = 'Confidence region'), linetype = 'dashed') +
    geom_line(aes(x = band$x, y = m_vals, color = 'True m')) +
    scale_color_manual(values = c('Estimated' = 'blue', 'Confidence region' = 'red', "True m" = 'green')) +
    ggtitle(paste0('t0=0, c=', c, 'd=', d, 'mp_type=', mp_type)) +
    ylab('m(t0, x) and SCR') +
    theme(legend.position = 'top') + 
    labs(title = paste0('t0=0, c=', c, 'd=', d, 'mp_type=', mp_type), 
         subtitle = paste0('m_case=', m_case, ', eps_case=', eps_case))

#   View this issue from fixed x0=0 perspective 
res_x <- auto.SCR(X, c=c, d=d, b_time=b_time, b_timese=b_timese,
                mp_type=mp_type, type='fixx', ops='CV', fix_num = 0, per=0.1)
df <- res_x[[1]]
extract_band_x <- function(df) {
  n_point <- length(df$t) / 3
  t_grid  <- df$t[1:n_point]
  idx_e   <- df$order == "esti"
  idx_u   <- df$order == "upper"
  idx_l   <- df$order == "lower"
  list(t = t_grid,
       esti  = df$y[idx_e],
       upper = df$y[idx_u],
       lower = df$y[idx_l])
}
band <- extract_band_x(df)

m_true <- model_cases[[m_case]]$m
m_vals <- m_true(t = band$t, x = rep(0, length(band$t)), delta = delta)
ggplot() +
    geom_line(aes(x = band$t, y = band$esti, color = 'Estimated')) +
    geom_line(aes(x = band$t, y = band$upper, color = 'Confidence region'), linetype = 'dashed') +
    geom_line(aes(x = band$t, y = band$lower, color = 'Confidence region'), linetype = 'dashed') +
    geom_line(aes(x = band$t, y = m_vals, color = 'True m')) +
    scale_color_manual(values = c('Estimated' = 'blue', 'Confidence region' = 'red', "True m" = 'green')) +
    ggtitle(paste0('t0=0, c=', c, 'd=', d, 'mp_type=', mp_type)) +
    ylab('m(t0, x) and SCR') +
    theme(legend.position = 'top') + 
    labs(title = paste0('t0=0, c=', c, 'd=', d, 'mp_type=', mp_type), 
         subtitle = paste0('m_case=', m_case, ', eps_case=', eps_case))
