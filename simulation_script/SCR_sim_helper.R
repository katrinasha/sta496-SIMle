library(SIMle)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(parallel)
library(snow)

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
  X_prev <- ifelse(m_case!=3, x0, 1)
  for (t in 1:n) {
    u <- t/n
    mu <- m_fun(u, X_prev)
    sig <- s_fun(u, X_prev)
    X[t] <- mu + sig * eps[t]
    X_prev <- X[t]
  }
  X
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

auto_SCR_coverage <- function(n_rep = 500, n = 100, t_grid = seq(0, 1, length.out = 101),
                              c = 10, d = 10, b_time = 'Legen', b_timese = 'Legen',
                              mp_type = 'algeb', type = 'fixt', ops = 'CV', per=0.2,
                              r=1, s=1, upper=10, delta=0, m_case, eps_case = 'tvAR2',
                              n_cores = NULL, debug = TRUE) {
  
  # Start timing
  total_start_time <- Sys.time()
  
  # Safe maximum for socket connections
  MAX_SAFE_CORES <- 60  # Conservative limit well below 128 socket connection limit
  
  # Determine number of cores if not specified, respecting socket limits
  if (is.null(n_cores)) {
    n_cores <- min(detectCores() - 1, n_rep, MAX_SAFE_CORES)
  } else {
    n_cores <- min(n_cores, MAX_SAFE_CORES)  # Cap requested cores to safe limit
  }
  
  # Calculate chunking if needed
  chunk_size <- n_cores
  n_chunks <- ceiling(n_rep / chunk_size)
  
  # Print diagnostic information
  if (debug) {
    cat("\nParallel Computing Setup Information:\n")
    cat("----------------------------------------\n")
    cat("Total available cores:", detectCores(), "\n")
    cat("Originally requested cores:", n_cores, "\n")
    cat("Safe maximum cores (socket limit):", MAX_SAFE_CORES, "\n")
    cat("Final core count:", n_cores, "\n")
    cat("Number of replications:", n_rep, "\n")
    cat("Chunk size:", chunk_size, "\n")
    cat("Number of chunks:", n_chunks, "\n")
    cat("System information:\n")
    cat("OS:", Sys.info()["sysname"], "\n")
    cat("Hostname:", Sys.info()["nodename"], "\n")
    cat("Current ulimit -n:", system("ulimit -n", intern=TRUE), "\n")
    cat("----------------------------------------\n\n")
  }
  
  # Note: cluster creation and exports are done per-chunk below to avoid
  # leaving stray socket connections open. We perform package/function
  # exports inside each chunk's cluster so workers have the needed objects.
  
  # Define the worker function with timing
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
  
  # Run parallel computations with chunking and progress monitoring
  comp_start_time <- Sys.time()
  all_results <- vector("list", n_rep)
  
  for (chunk in 1:n_chunks) {
    chunk_start <- (chunk-1) * chunk_size + 1
    chunk_end <- min(chunk * chunk_size, n_rep)
    chunk_indices <- chunk_start:chunk_end
    
    if (debug) {
      cat(sprintf("\nProcessing chunk %d/%d (replications %d-%d) using %d cores...\n", 
                 chunk, n_chunks, chunk_start, chunk_end, n_cores))
    }
    
  # Create fresh cluster for each chunk to avoid socket accumulation
  cl <- makeCluster(n_cores, type = "SOCK")
    
    # Export necessary objects to new cluster
    clusterExport(cl, c(
      "sim_x", "auto.SCR", "extract_band", "m1", "m2", "m3",
      "s1", "s2", "s3",
      "tvAR2_eps", "setar_eps", "bilinear_eps",
      "a1_fun", "a2_fun",
      "n", "delta", "m_case", "eps_case", "t_grid",
      "c", "d", "b_time", "b_timese", "mp_type", "type", "ops", "per"
    ), envir = environment())
    
    # Load required packages on new cluster
    clusterEvalQ(cl, {
      library(SIMle)
      library(data.table)
    })
    
    # Process this chunk
    chunk_results <- parLapply(cl, chunk_indices, worker_fun)
    all_results[chunk_indices] <- chunk_results
    
    # Clean up cluster after each chunk
    stopCluster(cl)
    
    if (debug) {
      cat(sprintf("Chunk %d completed. Time so far: %s\n", 
                 chunk, format(difftime(Sys.time(), comp_start_time))))
    }
  }
  
  # Extract results and timing information
  covered <- sapply(all_results, function(x) x$result)
  worker_times <- sapply(all_results, function(x) x$worker_time)
  
  if (debug) {
    cat("\nPerformance Statistics:\n")
    cat("----------------------------------------\n")
    cat("Total computation time:", difftime(Sys.time(), comp_start_time, units="mins"), "minutes\n")
    cat("Average time per replication:", mean(worker_times), "seconds\n")
    cat("Min time per replication:", min(worker_times), "seconds\n")
    cat("Max time per replication:", max(worker_times), "seconds\n")
    cat("Standard deviation:", sd(worker_times), "seconds\n\n")
    
    # Load distribution across nodes for all chunks
    node_counts <- table(sapply(all_results, function(x) x$node))
    cat("Load distribution across nodes:\n")
    print(node_counts)
    cat("\n")
  }
  
  # No global cluster remains; clusters were stopped after each chunk.
  
  # Print simulation parameters
  cat('\nSimulation Parameters:\n')
  cat('----------------------------------------\n')
  cat('Model:', m_case, '(', eps_case, ')\n')
  cat('Replications:', n_rep, '\n')
  cat('Sample size:', n, '\n')
  cat('Grid points:', length(t_grid), '\n')
  cat('Other parameters:', sprintf('c=%d, d=%d, type=%s\n', c, d, type))
  
  # Calculate and return coverage probability with confidence interval
  coverage_prob <- mean(covered)
  ci <- binom.test(sum(covered), n_rep)$conf.int
  
  list(
    coverage_probability = coverage_prob,
    confidence_interval = ci,
    total_time = difftime(Sys.time(), total_start_time, units="mins"),
    n_cores_used = n_cores,
    worker_times = worker_times
  )
  
}

# auto_SCR_coverage(n_rep=500, n=500, delta=1, m_case = 1, eps_case = 'tvAR2')
# auto_SCR_coverage(f_true=m1, delta=1, m_case = 1, eps_case = 'SETAR')
# auto_SCR_coverage(n_rep=1, n=10, delta=1, m_case = 3, eps_case = 'bilinear')
