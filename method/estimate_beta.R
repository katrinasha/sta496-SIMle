x <- read.csv("~/Documents/GitHub/sta496-SIMle/method/test_data.csv") %>% pull(x)
data <- data.frame(
    Y = x,
    t = seq(0, 1, length.out = length(x)), 
    X = c(NA, x[1:(length(x)-1)])
) %>%
    data.table()
data <- data[-1,] # first row has NA for x_{t-1}
m1 <- function(t, x, delta) (1 + x^2)^(-4) + delta * sin(2*pi*t) * exp(-x^2)

source('method/legendre.R')
# Estimate OLS regression coefficient in 
#     Y = \betaW + error
# where W is the design matrix constructed from basis functions

estimate_beta <- function(data, outcome, covariates, b_time, b_timese, c, d) {
    n <- nrow(data)
    Y <- data %>% pull(outcome)
    t <- data %>% pull(t)
    X <- data %>% pull(covariates)
    
    # construct design matrix W
    t_comp <- legendre_basis(c, t, type = 't')
    # ggplot(t_comp, aes(x = variable, y = value, color = as.factor(basis_order))) +
    #     geom_line() +
    #     labs(title = "Legendre basis functions for time", x = "t", y = "Basis function value", color = "Basis order")
    x_comp <- legendre_basis(d, X, type = 'x')

    W0 <- as.matrix(t_comp %>% pull(value) %>% matrix(nrow = n, ncol = c))
    W1 <- matrix(0, nrow = n, ncol = c * d)
    for (i in 1:n) {
        phi_t <- t_comp %>%
            filter(variable == t[i]) %>%
            arrange(basis_order) %>%
            pull(value)
        varphi_x <- x_comp %>%
            filter(variable == X[i]) %>%
            arrange(basis_order) %>%
            pull(value)
        W1[i, ] <- as.vector(t(outer(phi_t, varphi_x))) # phi_t1*varphi_x1, phi_t1*varphi_x2, ...
    }
    # W <- cbind(W0, W1[, -1]) # 1st column of W1 is all 1s, already included in W0; avoid collinearity
    W <- W1 # assumed no intercept in model to estimate

    # OLS estimate
    beta_hat <- solve(t(W) %*% W, tol = 1e-40) %*% t(W) %*% Y
    esti <- W %*% beta_hat
    return(list(beta_hat = beta_hat, esti = esti, W = cbind(W0, W1), w = x_comp))
}


#esti <- estimate_beta(data, outcome = 'Y', covariates = 'X', b_time = NULL, b_timese = NULL, c = 5, d = 5)
# true = m1(t = seq(0, 1, length.out = nrow(data)), x = x, delta = 1)
# data1 = data.frame(Y = true, t = data$t, X = x) %>% data.table()

true=m_true(t=rep(0, n_point), x=band$x, delta=1)
data1 = data.frame(Y = true[2:500], X = x[1:499], t=seq(0,1,length.out=499))

esti <- estimate_beta(data1, outcome = 'Y', covariates = 'X', b_time = NULL, b_timese = NULL, c = 15, d = 15)
df <- data.frame(esti = esti$esti, true = true[2:500], t = seq(0,1,length.out=499), x = data$X) %>%
    data.table()
ggplot(df, aes(x = true, y = esti)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    labs(title = "Estimated m(t,x) vs True m(t,x)", x = "True m(t,x)", y = "Estimated m(t,x)")
# +
#     geom_point(aes(y = true, color = "true")) +
#     geom_point(aes(y = esti, color = "esti")) +
#     labs(title = "Estimated m(t,x) vs True m(t,x)", x = "Time t", y = "m(t,x)", color = "Legend") +
#     theme(legend.position = "top") +
#     scale_color_manual(values = c("true" = "red", "esti" = "blue"), labels = c("True m(t,x)", "Estimated m(t,x)"))
