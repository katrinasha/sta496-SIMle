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
data <- data[-1,] # first row has NA for x_{t-1}

algeb_mapping_to_R <- function(y, s = 1) {
    # y in [-1, 1] -> R
    return(s * y / sqrt(1 - y^2))
}

algeb_mapping_to_unit <- function(x, s=1) {
    # x in R -> [-1, 1]
    return(x/sqrt(x^2 + s^2))
}

scr <- function(data, outcome, covariates, b_time, b_timese, c, d, m, fix_t, M = 1000, alpha = 0.05, C1 = 2000, C2 = 2000) {
    # c, d are set number of basis functions
    n <- nrow(data)

    # Step 1: Compute \hat{Pi}, estimate m(t,x) and residuals
    res <- estimate_beta(data, outcome = "Y", covariates = "X", b_time, b_timese, c, d)
    W <- res$W
    hatPi <- 1 / n * t(W) %*% W
    esti <- res$esti
    resid <- data$Y - esti

    # Setp 2: Construct a sequence of uniform grids on [0,1] x [0, 1]
    t_grid <- seq(0, 1, length.out = C1)
    y_grid <- seq(0, 1, length.out = C2)
    x_grid <- algeb_mapping_to_R(2 * y_grid - 1) # map to R

    basis_ti <- matrix(nrow = C1, ncol = c)
    for (l1 in (1:c)){
        basis_ti[, l1] <- legendre_basis(l1, t_grid, type = "t") %>%
            filter(basis_order == (l1-1)) %>%
            pull(value)
    }

    # Step 3: For each pair of (t,x) calculate sd --> \hat{h} (t,x)
    # \hatT(t,x) = \Xi \hatPi \hat bar_r \barL
    #   (i) compute \Xi (the boostrap)
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
    b <- matrix(nrow = c*d)
    xi <- seq(-10, 10, length.out = C1)
    t_idx <- unique(which(abs(fix_t - t_grid) == min(abs(fix_t - t_grid))))[[1]]
    basis_xi <- matrix(nrow=C1, ncol=c*d)
    for (i in 1:C1) {
        ti <- matrix(as.numeric(basis_ti[t_idx, ]), ncol = 1)
        xj <- legendre_basis(d, xi[i], type = "x")$value
        b <- cbind(b, kronecker(ti, xj))
        basis_xi[i,] <- xj
    }
    b <- b[, -1]
    
    hatT_list <- c()
    for (k in (1:M)) {
        hatT <- t(Xi_list[[k]][(c+1):(c+c*d)]) %*% solve(hatPi, tol = 1e-40) %*% diag(c * d) %*% b
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
    b1 <- matrix(nrow = C1, ncol = 1)
    for (l1 in 1:c) {
        for (l2 in 1:d) {
            b1 <- cbind(b1, basis_ti[t_idx, ] * matrix(basis_xi[, l2], ncol = 1))
        }
    }
    b1 <- b1[, -1]
    esti <- b1 %*% res$beta_hat

    return(data.frame(esti = esti, lower = esti - bound, upper = esti + bound))
}

#res <- scr(data, 'Y', 'X', b_time='Legen', b_timese='Legen', c=5, d=5, m=3, fix_t=0.5)
