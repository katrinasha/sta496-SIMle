library(dplyr)
library(data.table)
library(ggplot2)

source("method/estimate_beta.R")
x <- read.csv("~/Documents/GitHub/sta496-SIMle/method/test_data.csv") %>% pull(x)
data <- data.frame(
    Y = x,
    t = seq(0, 1, length.out = length(x)), 
    X = c(NA, x[1:(length(x)-1)])
) %>%
    data.table()
data <- data[-1, ] # first row has NA for x_{t-1}
m1 <- function(t, x, delta) (1 + x^2)^(-4) + delta * sin(2*pi*t) * exp(-x^2)

algeb_mapping_to_R <- function(y, s = 1) {
    # y in [-1, 1] -> R
    return(s * y / sqrt(1 - y^2))
}

algeb_mapping_to_unit <- function(x, s = 1) {
    # x in R -> [-1, 1]
    return(x / sqrt(x^2 + s^2))
}
vartheta_t <- function(x, t0, phi_x_list, K, h_t) {
  n <- length(x)
  stopifnot(j >= 1, j < n)
  
  # times on [0,1]
  t_grid <- (1:n) / n
  
  # indices where X_{i-j} exists
  idx <- 2:n
  
  # kernel weights in time around t0
  u <- (t_grid[idx] - t0) / h_t
  w_raw <- K(u)
  if (all(w_raw == 0)) {
    stop("All kernel weights are zero. Try a larger bandwidth h_t.")
  }
  w <- w_raw / sum(w_raw)
  
  d_j <- length(phi_x_list)
  vartheta_hat <- numeric(d_j)
  
  # compute E[phi_ell(X_{i-j}) | i/n near t0] via weighted average
  for (ell in 1:d_j) {
    phi_vals <- sapply(x[idx - j], phi_x_list[[ell]])   # phi_ell(X_{i-j})
    vartheta_hat[ell] <- sum(w * phi_vals)
  }
  
  vartheta_hat
}

scr <- function(data, outcome, covariates, b_time, b_timese, c, d, m, fix_t, M = 1000, alpha = 0.05, n_point = 2000) {
    # c, d are set number of basis functions
    n <- nrow(data)

    # Step 1: Compute \hat{Pi}, estimate m(t,x) and residuals
    res <- estimate_beta(data, outcome = "Y", covariates = "X", b_time, b_timese, c, d)
    W <- res$W
    hatPi <- 1 / n * t(W) %*% W
    esti <- res$esti
    resid <- data$Y - esti
    X <- data$X



    # Setp 2: Construct a sequence of uniform grids on [0,1] x [0, 1]
    t_grid <- seq(0, 1, length.out = n_point)
    t_idx <- unique(which(abs(fix_t - t_grid) == min(abs(fix_t - t_grid))))[[1]]
    y_grid <- seq(0, 1, length.out = n_point)
    # x_grid <- algeb_mapping_to_R(2 * y_grid - 1) # map to R
    x_grid <- seq(-10, 10, length.out = n_point)

    # Compute \hat \bar r * \bar I_1
    ## \hat \bar r
    ### compute b(t0,x)
    ti <- t_grid[t_idx]
    basis_ti <- matrix(nrow = c, ncol = 1)
    basis_ti <- legendre_basis(c, ti, type='t')$value # phi_1(ti),...,phi_c(ti)
    b_ti_x <- matrix(nrow = n_point, ncol = c * d)
    
    for (i in 1:length(x_grid)) {
        xi <- x_grid[i]
        basis_xi <- legendre_basis(d, xi / 10, type = NA)$value
        b_ti_xi <- kronecker(basis_ti, basis_xi)
        b_ti_x[i, ] <- b_ti_xi
    }
    f_ti <- matrix(nrow = n_point, ncol = c * d)
    var_theta_ti <- matrix(nrow = 1, ncol = d)
    # assumed lag r=1
    for (l2 in (1:d)) {
        # E phi_l2(X_1t)
        # use Gaussian kernel regression
        h <- n^(-1 / 5)
        idx <- 2:n
        u <- (seq(0, 1, length.out = n - 1) - fix_t) / h
        w <- dnorm(u) / sum(dnorm(u))
        phi_x <- legendre_basis(l2, X[idx - 1] / 10, type = NA) %>%
            filter(basis_order == l2 - 1) %>%
            pull(value)
        var_theta_ti[1, l2] <- sum(w * phi_x)
    }
    f_ti <- as.vector(t(kronecker(basis_ti, var_theta_ti)))

    r_ti_x <- b_ti_x - f_ti
    r_I <- cbind(matrix(0, ncol=c, nrow=n_point), r_ti_x)
    

 

    # Step 3: For each pair of (t,x) calculate sd --> \hat{h} (t,x)
    # \hatT(t,x) = \Xi \hatPi \hat bar_r \barL
    #   (i) compute \Xi (the boostrap)
    m=3
    eps_im_list <- list() # below (4.25)
    ui_list <- list() # (4.25)
    u_mi_list <- list() # below (4.25)
    for (i in 1:length(t_grid)) {
        ti <- t_grid[i]
        for (xl in x_grid[2:(length(x_grid) - 1)]) {
            w_ik <- legendre_basis(c, xl, type = "x") %>% pull(value)
            hat_u_ik <- w_ik * resid[i]
            ui_list[[i]] <- hat_u_ik
        }
    }
    # compute eps_im_list and u_mi_list
    for (i in 1:(length(t_grid) - m)) {
        hat_eps_im <- sum(resid[i:(i + m)])
        eps_im_list[[i]] <- hat_eps_im

        hat_u_mi <- 0
        for (o in c(i:(i + m))) {
            hat_u_mi <- hat_u_mi + ui_list[[o]]
        }
        u_mi_list[[i]] <- hat_u_mi
    }

    U_im_list <- list() # above (4.26)
    for (i in 1:(length(t_grid) - m)) {
        ti <- t_grid[i]
        U_im <- numeric(c + c * d)
        phi0_t <- legendre_basis(c, ti, type = "t") %>% pull(value)
        U_im[1:c] <- outer(eps_im_list[[i]], phi0_t)
        # phi1_t <-
        U_im[(c + 1):(c + c * d)] <- as.vector(t(outer(u_mi_list[[i]], phi0_t)))
        U_im_list[[i]] <- U_im
    }

    Xi_list <- list()
    # construct bootstrap samples
    for (b in 1:M) {
        set.seed(b)
        Ri_list <- rnorm(length(t_grid) - m, mean = 0, sd = 1)
        Xi <- 0
        for (i in 1:length(U_im_list)) {
            Xi <- Xi + U_im_list[[i]] * Ri_list[i]
        }
        Xi_list[[b]] <- 1 / sqrt((n - m - 1) * m) * Xi # (4.26): r=1

    }

    # Step 4: For each pair of (t,x), generate M many boostrap copies of \hat{T}(t,x)
    #     \bar r (4.5)
    # r = b-f
    
    hatT_list <- c()
    for (k in (1:M)) {
        hatT <- t(Xi_list[[k]][(c+1):(c+c*d)]) %*% solve(hatPi, tol = 1e-40) %*% r_I
        hatT_list[[k]] <- hatT
    }

    #     Compute order statistics \Tau_k
    df_T <- matrix(unlist(hatT_list), ncol = M)
    h_hat <- apply(df_T, 1, sd)

    tau_list <- c()
    for (k in 1:M) {
        tau_list[k] <- max(abs(df_T[, k]) / h_hat)
    }
    
    #     Compute critical value: the (1-alpha) quantile c(1-alpha)
    tau_ordered <- sort(tau_list)
    c_alpha <- tau_ordered[M*(1-alpha)]

    # Step 5: SCR m(t,x) +- c(1-alpha)\sqrt(n) \hat{h}(t,x)
    bound <- c_alpha / sqrt(n) * h_hat
    
    # Final: prepare result, estimate
    b1 <- matrix(nrow = n_point, ncol = 1)
    for (l1 in 1:c) {
        for (l2 in 1:d) {
            b1 <- cbind(b1, basis_ti[t_idx, ] * matrix(basis_xi[, l2], ncol = 1))
        }
    }
    b1 <- b1[, -1]
    esti <- b1 %*% res$beta_hat

    return(data.frame(x_grid = xi, esti = esti, lower = esti - bound, upper = esti + bound))
}

# res <- scr(data, 'Y', 'X', b_time='Legen', b_timese='Legen', c=5, d=5, m=3, fix_t=0.5)
true_fix_t <- m1(t = 0.5, x = res$x_grid, delta = 1)
# res %>%
#     mutate(true = true_fix_t) %>%
#     ggplot(aes(x = x_grid)) +
#     geom_line(aes(y = true, color = "true")) +
#     geom_line(aes(y=esti, color='esti'))
#     geom_line(aes(y = lower, color = "bound"), linetype='dashed') +
#     geom_line(aes(y = upper, color = "bound"), linetype='dashed') +
#     scale_color_manual(values=c('true'='black', 'bound'='red', 'esti'='blue'))
