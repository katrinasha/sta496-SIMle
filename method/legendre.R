# Legendre basis functions:
# 1. calculate coefficients in front of each order
# 2. compute the value of basis functions at given x

move_order <- function(poly){ # get the coefficient when polynomial multiple 1 times of x
  return(c(0,poly))
}

legendre_coef <- function(n){
    # n: number of Legendre basis functions
    # x: input value in [-1, 1] 
    coef_list <- list()
    coef_list[[1]] <- c(1)
    if (n ==1) {
        return(coef_list[[1]])
    } else if (n==2) {
        coef_list[[2]] <- c(0, 1)
        #return(coef_list[[2]]) # return (coef_list)
    } else {
        coef_list[[2]] <- c(0, 1)
        for (i in 3:n) {
            # np1 <- i+1
            # coef_list[[np1]] <-c(np1 + i)/np1 * move_order(coef_list[[i]]) - c(i/np1) * c(coef_list[[i-1]], 0, 0)

            coef_list[[i]] <- c((2*i)-3)/(i-1) * move_order(coef_list[[i-1]]) - c((i-2)/(i-1) * coef_list[[i-2]], 0, 0)
        }
        #return(coef_list[[n]]) # return(coef_list)
    }
    return(coef_list)
    

}
legendre_coef(5)

legendre_basis <- function(c, z, type = c('t','x')) {
    # c: number of basis functions (order = c-1)
    # x: input value in [-1,1]
    coef_list <- legendre_coef(c)

    # map x to [-1,1]
    if (type == 't') {
        x_mapped <- map_t_to_legendre(z)
    } else if (type == 'x') {
        x_mapped <- map_x_to_legendre(z)
    } else {
        stop("Invalid type. Must be 't' or 'x'.")
    }

    # compute each basis function at x_mapped
    basis_values <- list(c)
    for (i in 1:c) {
        poly_i <- coef_list[[i]]
        sum_i <- 0
        for (j in 1:length(poly_i)) {
            sum_i <- sum_i + poly_i[j] * (x_mapped)^(j-1)
        }
        basis_values[[i]] <- sum_i
    }
    n <- length(x_mapped)
    basis_values_df <- data.frame(
        basis_order = rep(0:(c-1), each = n),
        variable           = rep(z, times = c),
        value       = unlist(basis_values, use.names = FALSE)
        ) %>%
        data.table()
    
    return(basis_values_df)
}

map_t_to_legendre <- function(t) {
    # map t in [0,1] to x in [-1,1]
    return(2*t - 1)
}

map_x_to_legendre <- function(x, s=1) {
    # map x from R to [-1,1]
    x_temp <- x/(2*sqrt(x^2+s^2)) + 1/2 # now in [0,1]
    return(map_t_to_legendre(x_temp))
}

legendre_basis(3, c(0, 0.1, 1), type = 't')
legendre_basis(6, seq(0, 1, by=0.01), type = 't') %>%
    ggplot(aes(x = variable, y = value, color = as.factor(basis_order))) +
    geom_line() +
    labs(color = "Basis Order") +
    ggtitle("Legendre Basis Functions")
    
legendre_basis(6, seq(-3, 3, by=0.1), type = 'x') %>%
    ggplot(aes(x = variable, y = value, color = as.factor(basis_order))) +
    geom_line() +
    labs(color = "Basis Order") +
    ggtitle("Legendre Basis Functions Mapped from X in R")
