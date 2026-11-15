# data generative functions
m <- function(t, x, delta) {
    y_raw <- (x - 0.5)^2
    c <- as.integer(y_raw <= 2)
    return(c*y_raw)
    # 0
}

s <- function(t, x) {
    return(1)
}
# tvAR2
a1_fun <- function(u) 0.3
a2_fun <- function(u) 0.3 * sin(2 * pi * u)

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

case1 <- list(m = m, s = s)
model_cases <- list(case1)

sim_x <- function(n, delta = 1, m_case, eps_case, x0 = 0, burn = 500L) {
    stopifnot(m_case %in% 1:3)
    case_function <- model_cases[[m_case]]
    m_fun <- case_function$m
    s_fun <- case_function$s
    eps_fun <- errors[[eps_case]]
    print("m_fun=")
    print(m_fun)
    print("s_fun=")
    print(s_fun)

    # Generate x_t
    eps <- eps_fun(n, burn = burn)

    X <- numeric(n)
    x_prev <- x0
    data_used <- data.frame()
    for (t in 1:n) {
        u <- t / n
        m_val <- m_fun(u, x_prev, delta)
        s_var <- s_fun(u, x_prev)
        x_t <- m_val + s_var * eps[t]
        # x_t <- m_val
        data_used <- rbind(data_used, data.frame(t = u, x = x_t, x_prev = x_prev, eps = eps[t]))
        X[t] <- x_t
        x_prev <- x_t
    }
    data_used %>% data.table()
    return(X)
}

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

plot_fix_t_scr <- function(res, m_vals) {
df <- res[[1]]
band <- extract_band(df)
ggplot() +
    geom_line(aes(x = band$x, y = band$esti, color = "Estimated")) +
    geom_line(aes(x = band$x, y = band$upper, color = "Confidence region"), linetype = "dashed") +
    geom_line(aes(x = band$x, y = band$lower, color = "Confidence region"), linetype = "dashed") +
    geom_line(aes(x = band$x, y = m_vals, color = "True m")) +
    scale_color_manual(values = c("Estimated" = "blue", "Confidence region" = "blue", "True m" = "red")) +
    ggtitle(paste0("t0=", t, "c=", c, "d=", d, "mp_type=", mp_type)) +
    ylab("m(t0, x) and SCR") +
    theme(legend.position = "top") +
    labs(
        title = paste0("t0=", t, "c=", c, "d=", d, "upper=", upper, "mp_type=", mp_type),
        subtitle = paste0("m_case=", m_case, ", eps_case=", eps_case)
    ) +
    ylim(-5,5)

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
    t_grid <- seq(0, 1, length.out = 21)
    for (sim in 1:n_rep) {
        cat('Simulation ', sim, ' / ', n_rep, '\n')

        set.seed(sim)
        X <- sim_x(n, delta = delta, m_case = m_case, eps_case = eps_case, x0=0)
        ok_all_t <- TRUE

        for (t0 in t_grid) {
            cat('\tt0 = ', t0, '\n')

            res <- auto.SCR(X, c=c, d=d, b_time=b_time, b_timese=b_timese,
                            mp_type=mp_type, type='fixt', ops=ops, fix_num = t0, per=per, upper=upper)
            df <- res[[1]]
            band <- extract_band(df)
            m_true <- model_cases[[m_case]]$m
            m_vals <- m_true(t = rep(t0, length(band$x)), x = band$x, delta = delta)
            # plot(x=band$x, y=m_vals)
            # Plot SCR and true m at this t0
            p <- ggplot() +
                geom_line(aes(x = band$x, y = band$esti, color = "Estimated")) +
                geom_line(aes(x = band$x, y = band$upper, color = "Confidence region"), linetype = "dashed") +
                geom_line(aes(x = band$x, y = band$lower, color = "Confidence region"), linetype = "dashed") +
                geom_line(aes(x = band$x, y = m_vals, color = "True m")) +
                scale_color_manual(values = c("Estimated" = "blue", "Confidence region" = "red", "True m" = "green")) +
                ggtitle(paste0("Simulation ", sim, ", t0=", round(t0, 3), "c=", c, "d=", d, "mp_type=", mp_type)) +
                ylab("m(t0, x) and SCR") +
                theme(legend.position = "top") +
                labs(
                    title = paste0("Simulation ", sim, ", t0=", round(t0, 3), "c=", c, "d=", d, "mp_type=", mp_type),
                    subtitle = paste0("m_case=", m_case, ", eps_case=", eps_case)
                ) +
                ylim(-5,5)
            print(p)

            covered_t <- all(m_vals >= band$lower & m_vals <= band$upper)
            if (!covered_t) {
                ok_all_t <- FALSE
                break
            }
        }
        covered[sim] <- ok_all_t
        cat('Current uncovered vs. covered = ', table(covered), '\n')
    }
    return(mean(covered))

}
