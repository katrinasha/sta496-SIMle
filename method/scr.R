source("estimate_beta.R")

scr <- function(x, b_time, b_timese, c, d, B = 500, alpha = 0.05) {
    # c, d are set number of basis functions
    
    # Step 1: Compute \hat{Pi}, estimate m(t,x) and residuals

    # Setp 2: Construct a sequence of uniform grids on [0,1] x [0, 1]

    # Step 3: For each pair of (t,x) calculate sd --> \hat{h} (t,x)

    # Step 4: For each pair of (t,x), generate B many boostrap copies of \hat{T}(t,x)

    #     Compute order statistics \Tau_k
    #     Compute critical value: the (1-alpha) quantile c(1-alpha)

    # Step 5: SCR m(t,x) +- c(1-alpha)\sqrt(n) \hat{h}(t,x)
}