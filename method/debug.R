########## Step 1: use old version & simple m function ###########
library(Sie2nts)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(parallel)
library(snow)
# source('/Users/siyuansha/Documents/GitHub/sta496-SIMle/R/auto.SCR.R')
# source('/Users/siyuansha/Documents/GitHub/sta496-SIMle/R/SIMle.Legen.v1.R')
# source("/Users/siyuansha/Documents/GitHub/sta496-SIMle/method/SIMle.original_copy.R")
source('/Users/siyuansha/Documents/GitHub/sta496-SIMle/method/debug_helper.R')
# library(SIMle)

m_case = 1
eps_case='tvAR2'
n = 500
sim=1
set.seed(sim)
X <- sim_x(n = n, m_case = m_case, eps_case = eps_case)
eps <- tvAR2_eps(n = n)

c = 5
d = 5
b_time = "Legen"
b_timese = "Legen"
s = 1
r = 1
n_point = 2000
upper = 0.5
t = 0
mp_type = "algeb"
ops = "CV"
per=0.1
# res <- fix_t_scr(X,
#     t = t, c = c, d = d, m = "MV", b_time = b_time, b_timese = b_timese, mp_type = mp_type, r = r, s = s,
#     n_point = n_point, upper = upper
# )
res <- auto.SCR(X, c=c, d=d, b_time=b_time, b_timese=b_timese,
                            mp_type=mp_type, type='fixt', ops=ops, fix_num = t, per=per, upper=upper)
# m_vals <- model_cases[[m_case]]$m(t=rep(t, n_point), x=seq(-upper,upper,length.out=n_point))
# print(plot_fix_t_scr(res, m_vals = m_vals))
df <- res[[1]]
band <- extract_band(df)
m_true <- model_cases[[m_case]]$m
m_vals <- m_true(t = rep(t0, length(band$x)), x = band$x, delta = delta)

p <- ggplot() +
    # Truth (make bold + thicker + distinct color)
    geom_line(aes(x = band$x, y = m_vals, color = "Truth"), linewidth = 1.2) +

    # Estimation (medium thickness + dark color)
    geom_line(aes(x = band$x, y = band$esti, color = "Estimation"), linewidth = 1) +

    # Confidence band (lighter color + dashed + thinner)
    geom_line(aes(x = band$x, y = band$upper, color = "95% CI"), 
              linetype = "dashed", linewidth = 0.7) +
    geom_line(aes(x = band$x, y = band$lower, color = "95% CI"), 
              linetype = "dashed", linewidth = 0.7) +

    scale_color_manual(
        values = c(
            "Truth"       = "#D81B60",      # strong pink-red
            "Estimation"  = "#1F1F1F",      # deep charcoal black
            "95% CI"      = "#1E88E5"       # bright blue
        )
    ) +
    labs(
        x = "x",
        y = "f(x)",
        color = "",
        title = NULL
    ) +
    ylim(-1, 1) +
    theme_minimal(base_size = 13) +
    theme(
        legend.position = "top",
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    )
print(p)
ggsave('~/Downloads/ci_interval.png', p, width=6, height=3)
# plot(x=seq(-upper, upper, length.out=n_point), y=m(NA, seq(-upper, upper, length.out=n_point), NA))
# x_prev = 0
# x_list <- c(x_prev)
# n=100
# for (i in (1:n)){
#     x <- m(t = NA, x = x_prev, delta = NA)
#     x_list <- c(x_list, x)
#     x_prev=x
# }
# plot(x_list)
# fix fix_t_scr: no issue, same result as new implementation
# esti <- fix_t_esti(X, c, d, t, b_time, b_timese, mp_type, r, s, n_point, upper)
# length(esti[[1]])