library(SIMle)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)

# Setup
# model: Yi = m(ti, Xi)+sigma(ti,X_i-1)epsilon_i
# m(t,x) = (1+x^2)^{-4} + deltasin(2pit)exp(-x^2)

delta <- 1

m <- function(t, x, delta) {
  (1+x^2)^{-4} + delta*(2*pi*t)*exp(-x^2)
}

## ------------------------------------------------------------
## Helpers: time-varying coefficient functions a1(u), a2(u)
## (u = i/n in [0,1]). You can change these to match the paper.
a1_fun <- function(u) 0.3
a2_fun <- function(u) 0.3 * sin(2*pi*u)

## ------------------------------------------------------------
## Option A: epsilon_t is time-varying AR(2)
##   eps_t = a1(t)*eps_{t-1} + a2(t)*eps_{t-2} + eta_t
tvAR2_eps <- function(n, a1 = a1_fun, a2 = a2_fun, burn = 500L) {
  N <- n + burn
  eps <- rep(0, N)
  eta <- rnorm(N)
  for (t in 3:N) {
    u <- t/N
    eps[t] <- a1(u)*eps[t-1] + a2(u)*eps[t-2] + eta[t]
  }
  tail(eps, n)
}

## ------------------------------------------------------------
## Option B: epsilon_t is SETAR(1)
##   eps_t = a1(t)*eps_{t-1} + eta_t  if eps_{t-1} >= 0
##         = a2(t)*eps_{t-1} + eta_t  if eps_{t-1} <  0
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

## ------------------------------------------------------------
## Option C: epsilon_t is first-order bilinear
## A common specification (consistent with the text):
##   eps_t = (a1(t)*eta_t + a2(t)*eta_{t-1}) * eps_{t-1} + eta_t
bilinear_eps <- function(n, a1 = a1_fun, a2 = a2_fun, burn = 500L) {
  N <- n + burn
  eps <- rep(0, N)
  eta <- rnorm(N)
  for (t in 2:N) {
    u <- t/N
    eps[t] <- (a1(u)*eta[t] + a2(u)*eta[t-1]) * eps[t-1] + eta[t]
  }
  tail(eps, n)
}

## ------------------------------------------------------------
## m(t, x) and sigma(t, x) choices from the paper’s Section 5.1
## (Enter exactly as specified in the paper; below are templates.)
## t in [0,1], x in R, delta in [0,1].

m1 <- function(t, x, delta) (1 + x^2)^(-4) + delta * sin(2*pi*t) * exp(-x^2)
s1 <- function(t, x)        1.5 * exp(-x^2/2) * (2 + sin(2*pi*t))

m2 <- function(t, x, delta) (delta * (sin(2*pi*t) + 1)) * exp(-x^2/2)
s2 <- function(t, x)        0.5 * exp(-x^2) * (cos(2*pi*t) + 1)

## Piecewise sigma example shown in the text for setup (3):
m3 <- function(t, x, delta) 2*delta*(exp(-2*x^2) + abs(x)^(-1/2)*exp(-x^2/2))
s3 <- function(t, x) {
  ax <- abs(x)
  if (ax <= 1)     return(0.7 * (1 + x^2))
  if (ax <= 0.5)   return(1.4)            # <-- adjust breakpoints if needed
  return(2.0)
}

## ------------------------------------------------------------
## Main generator: X_t = m(t/n, X_{t-1}) + sigma(t/n, X_{t-1}) * eps_t
sim_x <- function(n, delta = 0, 
                  m_case = 1L, eps_case = "tvAR2",
                  a1 = a1_fun, a2 = a2_fun, burn = 500L, x0 = 0) {
  stopifnot(m_case %in% 1:3)
  # Choose m and sigma
  m_fun <- switch(m_case, 
                  function(t,x) m1(t,x,delta),
                  function(t,x) m2(t,x,delta),
                  function(t,x) m3(t,x,delta))
  s_fun <- switch(m_case, s1, s2, s3)
  
  # Generate epsilon_t
  eps <- switch(eps_case,
                "tvAR2"   = tvAR2_eps(n, a1, a2, burn),
                "SETAR"   = setar_eps(n, a1, a2, burn),
                "bilinear"= bilinear_eps(n, a1, a2, burn),
                stop("eps_case must be one of: 'tvAR2','SETAR','bilinear'"))
  
  # Generate X_t
  X <- numeric(n)
  X_prev <- x0
  for (t in 1:n) {
    u <- t/n
    mu <- m_fun(u, X_prev)
    sig <- s_fun(u, X_prev)
    X[t] <- mu + sig * eps[t]
    X_prev <- X[t]
  }
  X
}

n <- 100 #500
n_sim <- 5000
ts <- sim_x(n)
ts_df <- data.frame(xt = ts, x_tm1 = lag(ts), t = seq(0, 1, length.out = n))
ggplot(ts_df, aes(x=x_tm1, y = xt)) + geom_line()
# plot of the TS
plot(ts, type = 'l')

# run method
res <- auto.SCR(ts, c=10, d=10, b_time = 'Legen', b_timese = 'Legen', mp_type='algeb', type = 'fixt', ops = 'CV', fix_num = 0, per = 0.9)
str(res) # 1 data frame
head(res[[1]])


# put result into long format
df <- data.frame(esti = res[[1]] %>% filter(order == 'esti') %>% pull(y), 
           upper = res[[1]] %>% filter(order == 'upper') %>% pull(y),
           lower = res[[1]] %>% filter(order == 'lower') %>% pull(y)); head(df)
dim(df) # Q: what does each row represent?

# plot the result
df %>%
  mutate(row_num = seq_len(n())) %>%
  ggplot(aes(x = row_num)) +
  geom_line(aes(y = esti), color = "blue") +
  geom_line(aes(y = lower), color = "red", linetype = "dashed") +
  geom_line(aes(y = upper), color = 'red', linetype = 'dashed') # Q: but this doesn't approximate the TS

reform_df <- function(res) {
  data.frame(esti = res[[1]] %>% filter(order == 'esti') %>% pull(y), 
             upper = res[[1]] %>% filter(order == 'upper') %>% pull(y),
             lower = res[[1]] %>% filter(order == 'lower') %>% pull(y))
}

plot_df <- function(df) {
  df %>%
    mutate(row_num = seq_len(n())) %>%
    ggplot(aes(x = row_num)) +
    geom_line(aes(y = esti), color = "blue") +
    geom_line(aes(y = lower), color = "red", linetype = "dashed") +
    geom_line(aes(y = upper), color = 'red', linetype = 'dashed') 
}

###### plot ground truth function #########
f1 <- function(t, x, delta=1) {
  m1(t, x, delta) + s1(t,x)
}
t <- seq(0, 1, length.out = 100)
x <- seq(-10, 10, length.out = 2000)
z <- outer(t, x, f1)
plot_ly(x = x, y = t, z = z) %>%
  add_surface()

ts <- sim_x(n)

x_tm1 <- c(NA, ts[1:n-1])
z <- outer(t, x_tm1, f1)
plot_ly(x = x_tm1, y = t, z = log(z)) %>%
  add_surface()

##### compute coverage ##########

fix_nums <- seq(0, 1, length.out = 10)

for (fix_num in fix_nums) {
  ts <- sim_x(n)
  scr <- auto.SCR(ts, c=10, d=10, b_time = 'Legen', b_timese = 'Legen', mp_type='algeb', type = 'fixt', ops = 'CV', fix_num, per = 0.9)
  # Q: what does output of scr represent in each coord
  plot_df(reform_df(scr))
}
