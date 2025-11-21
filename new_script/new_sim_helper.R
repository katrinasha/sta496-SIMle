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
  tail(eps, n)
#   eps[1:n]
    # eps
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
#   eps[1:n]
    # eps
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
#   eps[1:n]
    # eps
}

errors <- list(
  tvAR2 = tvAR2_eps,
  SETAR = setar_eps,
  bilinear = bilinear_eps
)

################### Main functions: m and sigma #####################
m1 <- function(t, x, delta) (1 + x^2)^(-4) + delta * sin(2*pi*t) * exp(-x^2)
s1 <- function(t, x)        1.5 * exp(-x^2/2) * (2 + sin(2*pi*t))
s_flat <- function(t, x) 1   # replace s1,s2,s3 temporarily
# case1 <- list(m = m1, s = s1)
case1 <- list(m=m1, s=s_flat)

m2 <- function(t, x, delta) (delta * sin(2*pi*t) + 1) * exp(-x^2/2)
s2 <- function(t, x)        0.5 * exp(-x^2) * cos(2*pi*t) + 1
case2 <- list(m = m2, s = s_flat)


m3 <- function(t, x, delta) 2*t*(delta*exp(-2*t*x^2)) + pi^{-1/2}*exp(-x^2/2)
s3 <- function(t, x) {
  ax <- abs(x)
  if (ax <= 1)     return(0.7 * (1 + x^2))
  if (ax > 1 & 0 <= t & t < 0.5)   return(1.4) 
  return(2.0)
}
case3 <- list(m = m3, s = s_flat)

case1 <- list(m = m1, s = s1)
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
# plot(sim_x(n, m_case = m_case, eps_case = eps_case), type='l')

################### Verify SCR coverage #####################
extract_band <- function(df) {
  n_point <- length(df$t) / 3
  x_grid  <- df$t[1:n_point]
  idx_e   <- df$order == "esti"
  idx_u   <- df$order == "upper"
  idx_l   <- df$order == "lower"
  list(x = x_grid,
       esti  = df$y[idx_e],
       upper = df$y[idx_u],
       lower = df$y[idx_l])
}

auto_SCR_coverage <- function(n_rep = 1000, n = 500, c = 5, d = 5, delta=1,
                            b_time = 'Legen', b_timese = 'Legen', mp_type = 'algeb', type = 'fixt', 
                            ops = 'CV', per = 0.2, m_case, eps_case) {
    covered <- logical(n_rep)
    t_grid <- seq(0, 1, length.out = 101)

    for (s in 1:n_rep) {
        cat('Simulation ', s, ' / ', n_rep, '\n')

        set.seed(s)
        X <- sim_x(n, delta = delta, m_case = m_case, eps_case = eps_case, x0=-10)
        ok_all_t <- TRUE
        for (t0 in t_grid) {
            cat('\t t0 = ', t0, '\n')
            # res <- auto.SCR(X, c = c, d = d, b_time = b_time, b_timese = b_timese,
            #                 mp_type = mp_type, type = type, ops = ops, fix_num = t0, per = per)
            res <- auto.SCR(X, c=c, d=d, b_time=b_time, b_timese=b_timese,
                            mp_type=mp_type, type='fixx', ops=ops, fix_num = t0, per=per)
            df <- res[[1]]
            band <- extract_band(df)
            m_true <- model_cases[[m_case]]$m
            m_vals <- m_true(t = rep(t0, length(band$x)), x = band$x, delta = delta)
            m_vals <- m_true(t = band$x, x = rep(t0, length(band$x)), delta = delta)

            ggplot() +
                geom_line(aes(x = band$x, y = band$esti), color = 'blue') +
                geom_line(aes(x = band$x, y = band$upper), linetype = 'dashed', color = 'red') +
                geom_line(aes(x = band$x, y = band$lower), linetype = 'dashed', color = 'red') +
                geom_line(aes(x = band$x, y = m_vals), color = 'green') +
                ggtitle(paste0('Simulation ', s, ', t0=', round(t0, 3), 'c=', c, 'd=', d, 'mp_type=', mp_type)) +
                ylab('m(t0, x) and SCR') +
                # xlab('x') +
                scale_color_manual(values = c('blue' = 'Estimated m', 'red' = 'SCR Bands', 'green' = 'True m')) +
                theme(legend.position = 'top') + 
                labs(title = paste0('Simulation ', s, ', t0=', round(t0, 3), 'c=', c, 'd=', d, 'mp_type=', mp_type), 
                     subtitle = paste0('m_case=', m_case, ', eps_case=', eps_case))
            #  t0 = t0 + 0.01

            covered_t <- all(m_vals >= band$lower & m_vals <= band$upper)
            if (!covered_t) {
                ok_all_t <- FALSE
                break
            }
        }
        covered[s] <- ok_all_t
        cat('Current coverage: FALSE vs. True = ', table(covered), '\n')

    }
    return(mean(covered))

}
t0=0
s=4
n_rep=10
n = 500
c=5
d=c
delta=1
b_time = 'Legen'
b_timese = 'Legen'
mp_type = 'algeb'
type = 'fixx'
ops = 'CV'
per = 0.2
m_case = 1
eps_case = 'tvAR2'

#auto_SCR_coverage(m_case = 1, eps_case = 'tvAR2', n_rep=10)

# plot a sample X in 3D
library(plotly)
X <- sim_x(n, delta = 1, m_case = 1, eps_case = 'tvAR2', x0=0)
df <- data.frame(
    t = seq(1, n)/n,
    x = X,
    x_prev = c(0, X[1:499])
)
plot_ly(df, x = ~t, y = ~x_prev, z = ~x, type = 'scatter3d', mode = 'points') %>%
    layout(scene = list(xaxis = list(title = 't'),
                        yaxis = list(title = 'X_{t-1}'),
                        zaxis = list(title = 'X_t')))
# add a surface for m(t, x_prev)
t_seq <- seq(0, 1, length.out = 50)
x_seq <- seq(-10, 10, length.out = 50)
m_vals <- matrix(0, nrow = length(t_seq), ncol = length(x_seq))
for (i in 1:length(t_seq)) {
    for (j in 1:length(x_seq)) {
        m_vals[i, j] <- m1(t_seq[i], x_seq[j], delta = 1)
    }

}
plot_ly(x = t_seq, y = x_seq, z = m_vals) %>%
    add_surface() %>%
    layout(scene = list(xaxis = list(title = 't'),
                        yaxis = list(title = 'X_{t-1}'),
                        zaxis = list(title = 'm(t, X_{t-1})')))

# X2 vs. X3
X2 = X


covered <- logical(n_rep)
t_grid <- seq(0, 1, length.out = 101)

n_rep=10
n=500
m_case=1
c=5
d=5
mp_type='algeb'
b_time = 'Legen'
b_timese = 'Legen'

result <- list()
for (s in 1:n_rep) {
    cat('Simulation ', s, ' / ', n_rep, '\n')

    set.seed(s)
    X <- sim_x(n, delta = delta, m_case = m_case, eps_case = eps_case, x0=0) %>% scale()
    #scr_x_0 <- auto.SCR(X, c=c, d=d, b_time=b_time, b_timese=b_timese,
                          mp_type=mp_type, type='fixx', ops=ops, fix_num = 0, per=per, s=1)
    scr_x_0 <- fix.SCR(X, c=c, d=d, m='MV', b_time=b_time, b_timese=b_timese,
                         mp_type=mp_type, type='fixx', fix_num = 0, r=1, s=1, n_point = 1000, upper=10)
    df <- scr_x_0[[1]]
    band <- extract_band(df)
    m_true <- model_cases[[m_case]]$m
    model_cases[[m_case]]$s
    m_vals <- m_true(x = rep(0, length(band$x)), t = band$x, delta = delta)
    p <- ggplot() +
        geom_line(aes(x = band$x, y = band$esti), color = 'blue') +
        geom_line(aes(x = band$x, y = band$upper), linetype = 'dashed', color = 'red') +
        geom_line(aes(x = band$x, y = band$lower), linetype = 'dashed', color = 'red') +
        geom_line(aes(x = band$x, y = scale(m_vals)), color = 'black') +
        ggtitle(paste0('Simulation ', s, ', x=0,c=', c, 'd=', d, 'mp_type=', mp_type)) +
        ylab('m(t, 0) and SCR') +
        # xlab('x') +
        #scale_color_manual(values = c('blue' = 'Estimated m', 'red' = 'SCR Bands', 'green' = 'True m')) +
        theme(legend.position = 'top') + 
        labs(title = paste0('Simulation ', s, ', x=0,c=', c, 'd=', d, 'b_time=', b_time, 'b_timese=', b_timese, 'mp_type=', mp_type), 
             subtitle = paste0('n=', n, ', m_case=', m_case, ', eps_case=', eps_case))
    p
    #ggplotly(p)
    covered <- all(m_vals >= band$lower & m_vals <= band$upper)
    sum(!(m_vals >= band$lower & m_vals <= band$upper))
    if (!covered) {
      print('Not covered')
      result[[paste0('s', s)]] <- p
    }

}

result
m_case=2
m1 <- function(t, x, delta) 0
case1 <- list(m = m1, s = s1)
model_cases <- list(case1, case2, case3)
