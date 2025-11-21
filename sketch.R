  worker_fun <- function(s) {
    worker_start_time <- Sys.time()
    set.seed(s)
    ok_all_t <- TRUE
    ts <- sim_x(n, delta, m_case = m_case, eps_case = eps_case)
    m_options <- list(m1, m2, m3)

    n_t0 <- length(t_grid[2:(length(t_grid)-1)])
    t0_times <- numeric(n_t0)
    t0_idx <- 1

    for (t0 in t_grid[2:(length(t_grid)-1)]) {
      t0_start <- Sys.time()
      res <- auto.SCR(ts, c=c, d=d, b_time=b_time, b_timese=b_timese,
                      mp_type=mp_type, type=type, ops=ops, t0, per=per)
      df <- res[[1]]
      band <- extract_band(df)
      f_true <- m_options[[m_case]]
      m_true <- f_true(t=rep(t0, length(band$x)), x=band$x, delta=delta)
      covered_t <- all(m_true >= band$lower & m_true <= band$upper)
      t0_times[t0_idx] <- as.numeric(difftime(Sys.time(), t0_start, units="secs"))
      t0_idx <- t0_idx + 1

      if (!covered_t) {
        ok_all_t <- FALSE
        break
      }
    }

    # Return both result and timing information
    list(
      result = ok_all_t,
      worker_time = as.numeric(difftime(Sys.time(), worker_start_time, units="secs")),
      t0_times = t0_times,
      node = Sys.info()["nodename"],
      pid = Sys.getpid()
    )
  }


  # Settings
B <- 1000               # number of Monte Carlo repetitions
n <- 500                # sample size for each simulated time-series
seed_start <- 123
type <- "fixt"          # or "fixx", depending on what fix.SCR expects
ops <- "CV"             # selection method: "AIC","BIC","CV","Kfold"
per <- 0.1              # if CV used
k <- 5                  # if Kfold used
r <- 1
s <- 1
b_time <- "Legen"         # example basis for time
b_timese <- "Legen"     # example basis for x
m_window <- "MV"        # SCR window selection
upper <- 10
fixed_t <- 0.5          # when type == "fixt", the t at which you test the SCR
x_grid <- seq(-upper, upper, length.out = 200)  # grid of x for checking coverage

# Functions the user must provide:
# model_sim(n) -> vector ts (length n): simulates a time series under the known truth
# evaluate_true(t_or_x_grid) -> numeric vector of true m(t,x) values (same ordering as x_grid)
# (If your truth is m(t,x) you will supply the correct slice: evaluate_true(fixed_t, x_grid))

# Storage
sim_cover_indicator <- logical(B)        # TRUE if SCR covers entire true curve on x_grid
sim_pointwise_count <- numeric(length(x_grid)) # counts of times each x grid point is covered
selected_cd <- matrix(NA, nrow = B, ncol = 2)
band_widths <- numeric(B)                # average band width (over grid) per sim

set.seed(seed_start)
for(b in seq_len(B)){
  ts_b <- model_sim(n)    # user-supplied simulator returning time series vector

  # Run auto.SCR (auto picks c,d then calls fix.SCR)
  # auto.SCR returns a list of dataframes; each dataframe has columns:
  #   first col = input values (x or t),
  #   second col = concatenated values (esti, upper, lower),
  #   third col = order factor (e.g. "esti","upper","lower")
  # For type = "fixt", fix.SCR likely returns a single data frame per lag/variates. 
  scr_res <- auto.SCR(ts_b, c = <c_max>, d = <d_max>,
                      b_time = b_time, b_timese = b_timese,
                      mp_type = "algeb", type = type, ops = ops,
                      m = m_window, fix_num = fixed_t, r = r, s = s,
                      per = per, k = k, upper = upper)

  # Optionally capture selected c,d (best_cd.auto.fit returns them before calling fix.SCR)
  # best_cd.auto.fit returns two values (c,d) and auto.SCR used them internally.
  # If auto.SCR does not return them, you can call best_cd.auto.fit(ts_b, ...) directly:
  best_cd <- best_cd.auto.fit(ts_b, c = <c_max>, d = <d_max>,
                              b_time = b_time, b_timese = b_timese,
                              mp_type = "algeb", ops = ops, per = per, k = k, r = r, s = s)
  selected_cd[b, ] <- as.numeric(best_cd[1:2])

  # Extract the SCR dataframe for the test (for type == "fixt" auto.SCR likely returns one df)
  # Suppose scr_res is a list of length r (one per lag) and we evaluate the first element:
  df <- scr_res[[1]]  # adapt if auto.SCR returns different structure

  # The returned df structure in this package (see functions fix_t_scr / fix_x_scr) uses:
  #   df$x (or df$t) and df$y with order factor repeated blocks ("esti","upper","lower")
  # If it's the same format:
  # split blocks by order
  esti_block  <- df[df$order == "esti", ]
  upper_block <- df[df$order == "upper", ]
  lower_block <- df[df$order == "lower", ]

  # ensure the x-grid matches; if not, interpolate
  x_out <- esti_block$x  # input values returned by auto.SCR
  est_vals <- esti_block$y
  upper_vals <- upper_block$y
  lower_vals <- lower_block$y

  # If x_out != x_grid, evaluate true on x_out (preferred)
  true_vals <- evaluate_true(fixed_t, x_out)  # user function: true m(t, x_out)

  # Check pointwise coverage
  covered_pointwise <- (true_vals >= lower_vals) & (true_vals <= upper_vals)
  sim_pointwise_count <- sim_pointwise_count + as.numeric(covered_pointwise)

  # Simultaneous coverage for the slice: true curve must be inside band for ALL x_out points
  sim_cover_indicator[b] <- all(covered_pointwise)

  # average band width (mean(upper - lower)) for diagnostic
  band_widths[b] <- mean(upper_vals - lower_vals)
} # end B loop

# Summaries
sim_simultaneous_coverage <- mean(sim_cover_indicator)            # empirical simultaneous coverage
sim_simultaneous_se <- sqrt(sim_simultaneous_coverage*(1-sim_simultaneous_coverage)/B)
pointwise_coverage <- sim_pointwise_count / B                    # vector over x_grid
avg_band_width <- mean(band_widths)
cd_selected_table <- as.data.frame(table(selected_cd[,1], selected_cd[,2])) # frequency of selected (c,d)

# Output
list(
  B = B,
  simultaneous_coverage = sim_simultaneous_coverage,
  simultaneous_se = sim_simultaneous_se,
  pointwise_coverage = data.frame(x = x_grid, coverage = pointwise_coverage),
  avg_band_width = avg_band_width,
  selected_cd_freq = cd_selected_table
)