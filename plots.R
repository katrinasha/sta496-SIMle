library(dplyr)
library(ggplot2)
library(data.table)
data <- read.csv("/Users/siyuansha/Documents/GitHub/sta496-SIMle/hkepd_dat.csv") %>%data.table()
head(data)

df <- data %>%
    mutate(date = seq.Date(from = as.Date("1994-01-01"),
         to   = as.Date("1995-12-31"),
         by   = "day")) %>%
    data.table() %>%
    melt(measure.vars = c("SO2", "NO2", "RSP", "Total", 'Temp', 'Circulatory'))
library(ggplot2)
library(patchwork)

# Shared theme for consistency
base_theme <- theme_minimal(base_size = 13) +
    theme(
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = element_text(size = 11, color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )

# --- Hospital Admissions Plot ---
p_admission <- df %>% 
    filter(variable == "Circulatory") %>% 
    ggplot(aes(x = date, y = value)) +
    geom_area(fill = "#377eb8", alpha = 0.25) +
    geom_line(color = "#377eb8", size = 1) +
    # geom_smooth(method = "loess", se = FALSE, color = "black", linetype = "dashed", size = 0.8) +
    labs(
        y = "Hospital admissions",
        x = "",
        title = "Daily Hospital Admissions for Circulatory Disease in Hong Kong",
        # subtitle = "Daily series with smoothed trend"
    ) +
    base_theme

# --- Dust Level Plot ---
p_RSP <- df %>% 
    filter(variable == "RSP") %>% 
    ggplot(aes(x = date, y = value)) +
    geom_area(fill = "#e41a1c", alpha = 0.25) +
    geom_line(color = "#e41a1c", size = 1) +
    # geom_smooth(method = "loess", se = FALSE, color = "black", linetype = "dashed", size = 0.8) +
    labs(
        y = "Dust level (RSP)",
        x = "Date",
        title = "Daily Dust Level in Hong Kong",
        # subtitle = "High particulate matter levels may influence hospital admissions"
    ) +
    base_theme

# Combine vertically
p_admission_RSP <- p_admission / p_RSP
ggsave('~/Downloads/HK_admission_RSP.png', p_admission_RSP, width=8, height=4.5)

p1 <- ggplot(df %>% filter(variable %in% c('RSP', "Total")), aes(x = date, y = value, fill = variable)) +
    geom_line() +
    facet_wrap(~variable, ncol=1); p1

    ggplot(df %>% filter(variable %in% c('RSP', "Total")), aes(x = RSP, y = Total, fill = variable)) +
    geom_point()


    facet_wrap(~variable, ncol=1)
ggplot(data, aes(x = RSP, y = Total)) +
    geom_line()

library(plotly)
# plot_ly(
#   data = data,              # replace with your data.frame/tibble
#   x = ~X,
#   y = ~RSP,
#   z = ~Total,
#   type = "scatter3d",
#   mode = "markers",
#   marker = list(size = 3, color = ~Total, colorscale = "Viridis", showscale = TRUE)
# ) |>
#   layout(scene = list(
#     xaxis = list(title = "X"),
#     yaxis = list(title = "RSP"),
#     zaxis = list(title = "Total")
#   ))

library(ggplot2)
library(mgcv)  # for GAM if not already loaded

# Shared base theme (same as before)
base_theme <- theme_minimal(base_size = 13) +
  theme(
    axis.title.x  = element_text(face = "bold", size = 12),
    axis.title.y  = element_text(face = "bold", size = 12),
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    panel.grid.minor = element_blank()
  )
data <- data %>% mutate(date = seq.Date(
    from = as.Date("1994-01-01"),
    to = as.Date("1995-12-31"),
    by = "day"
), month = month(date), season = case_when(
    month %in% c(12, 1, 2) ~ "Winter",
    month %in% c(3, 4, 5) ~ "Spring",
    month %in% c(6, 7, 8) ~ "Summer",
    month %in% c(9, 10, 11) ~ "Fall"
))


p_summer <- ggplot(data %>% filter(season %in% c("Summer", "Spring"), RSP <= 120, aes(x = RSP, y = Circulatory))) +
    # scatter points
    geom_point(alpha = 0.2, color = "grey40") +

    # nonlinear smooth (GAM) - blue, solid
    geom_smooth(
        method = "loess",
        # formula = y ~ s(x, k = 20),
        se = TRUE,
        color = "#377eb8", # same blue as earlier
        linewidth = 1
    ) +

    # linear fit - red, dashed
    geom_smooth(
        method    = "lm",
        se        = FALSE,
        color     = "#e41a1c", # same red as earlier
        linetype  = "dashed",
        linewidth = 0.8
    ) +
    labs(
        x = "Dust level (RSP)",
        y = "Hospital admissions",
        title = "Circulatory Disease Hospital Admissions vs. Dust Level",
        # subtitle = "Blue: flexible GAM smooth; Red dashed: linear fit"
    ) +
    xlim(c(10, 100)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(50, 170)) +
    base_theme
p_summer

  p_winter <- ggplot(data %>% filter(RSP <= 120, season %in% c("Winter", "Fall")), aes(x = RSP, y = Circulatory)) +
      # scatter points
      geom_point(alpha = 0.2, color = "grey40") +

      # nonlinear smooth (GAM) - blue, solid
      geom_smooth(
          method = "loess",
          #   formula = y ~ s(x, k = 10),
          se = TRUE,
          color = "#377eb8", # same blue as earlier
          linewidth = 1
      ) +

      # linear fit - red, dashed
      geom_smooth(
          method    = "lm",
          se        = FALSE,
          color     = "#e41a1c", # same red as earlier
          linetype  = "dashed",
          linewidth = 0.8
      ) +
      labs(
          x = "Dust level (RSP)",
          y = "Hospital admissions",
          title = "Circulatory Disease Hospital Admissions vs. Dust Level",
          # subtitle = "Blue: flexible GAM smooth; Red dashed: linear fit"
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(c(10, 100)) +
      ylim(c(50,170))+
      base_theme; p_winter
  

p_summer+p_winter
    ggsave("~/Downloads/winter_trend.png", p_winter, width = 8, height = 4)
  ggsave('~/Downloads/spring_trend.png', p_summer, width=8, height=4)
  
p2 <- ggplot(data, aes(x = RSP, y = Total)) +
  # scatter points
  geom_point(alpha = 0.3, color = "grey40") +
  
  # nonlinear smooth (GAM) - blue, solid
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    color   = "#377eb8", # same blue as earlier
    linewidth = 1
  ) +
  
  # linear fit - red, dashed
  geom_smooth(
    method    = "lm",
    se        = FALSE,
    color     = "#e41a1c",  # same red as earlier
    linetype  = "dashed",
    linewidth = 0.9
  ) +
  
  labs(
    x = "Dust",
    y = "Admission",
    title    = "Relationship between Daily Temperature and Dust Level",
    # subtitle = "Blue: flexible GAM smooth; Red dashed: linear fit"
  ) +
  base_theme; p2

ggsave('~/Downloads/HK_nonlinear.png', p2, width=8, height=3)
  p3 <- ggplot(data, aes(x = X, y = Total)) +
      geom_line() +
      ggtitle("non-stationary")
  p3
  
  ggplot(data, aes(x = X, y = Total)) +
      geom_line() +
      ggtitle("time-varying")
  
  # interaction between Total ~ RSP + time
  plot_ly(
  data = data,
  x = ~X,        # time
  y = ~RSP,             # air quality
  z = ~Total,      # hospital admissions
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Time"),
      yaxis = list(title = "RSP"),
      zaxis = list(title = "Admissions")
    ),
    title = "3D Relationship: Time, RSP, and Hospital Admissions"
  )



stock <- read.csv("~/Downloads/replication-files-garch-0-nonstationary/bloomberg/NYA-long Index.csv") %>% data.table()

stock <- read.csv("~/Downloads/replication-files-garch-0-nonstationary/bloomberg/data1/CPS US Equity.csv") %>% data.table()
ggplot(stock, aes(x=as.Date(date), y=PX_LAST))+geom_line()
#   data$prev_Total <- c(NA, data$Total[0:729])
#   auto_RSP
#   data <- data %>% mutate(all_part = SO2+NO2+RSP)

#   ggplot(data, aes(x = prev_RSP, y = RSP)) +
#     geom_point(alpha = 0.3) +
#     geom_smooth(formula = y ~ s(x, k = 30), method = "gam", se = FALSE, color = "green", span = 0.8) +
#     geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
#     # geom_smooth(method = 'loess', span=0.5, se=FALSE, color='blue')
#   labs(
#     x = "PM2.5 concentration",
#     y = "Respiratory admissions",
#     title = "Nonlinear relationship between PM2.5 and hospital admissions",
#     subtitle = "Blue: flexible smooth; Red dashed: linear fit"
#   ) +
#   theme_minimal()

# mod <- lm(Total ~ SO2 + NO2, data=data)
# resid(mod)
# data$resid = resid(mod)
# plot(data$resid)

########### basis functions #######
# Install if needed:
# install.packages(c("ggplot2", "splines", "tidyr", "dplyr"))

library(ggplot2)
library(splines)
library(tidyr)
library(dplyr)

# Grid on [0, 1]
x <- seq(0, 1, length.out = 200)

# Helper: build a B-spline basis and return a long data frame for plotting
make_bs_df <- function(x, df, degree = 3) {
  B <- bs(x, df = df, degree = degree, intercept = TRUE)
  colnames(B) <- paste0("B", seq_len(ncol(B)))
  
  as.data.frame(B) |>
    mutate(x = x) |>
    pivot_longer(
      cols = starts_with("B"),
      names_to = "basis",
      values_to = "value"
    ) |>
    mutate(df = df,
           degree = degree)
}

# Choose a few sieve sizes (numbers of basis functions)
df_list <- c(5)

bs_df <- bind_rows(lapply(df_list, function(df) make_bs_df(x, df = df, degree = 3)))

bs_df$df <- factor(bs_df$df, levels = df_list,
                   labels = paste(df_list, "basis functions"))

# Plot: B-spline basis functions for different sieve sizes
# p_bspline <- ggplot(bs_df, aes(x = x, y = value, color = basis)) +
#     geom_line(size = 1) +
#     facet_wrap(~df, ncol = 1) +
#     labs(
#         x = "Input",
#         y = "Basis function value"
#     ) +
#     theme_minimal() +
#     theme(legend.position = "none")
# p_bspline

df_legen <- legendre_basis(5, seq(-1, 1, by = 0.01), type = NA) %>%
    mutate(basis_order = as.factor(basis_order)) 
p_legen <- ggplot(df_legen, aes(x = variable, y = value, color = basis_order)) +
    geom_line() +
    labs(x='Input', y='Basis function value')+
    theme(legend.position = 'none')

# Install if needed:
# install.packages(c("ggplot2", "tidyr", "dplyr"))

library(ggplot2)
library(tidyr)
library(dplyr)

# Grid on [0, 1]
x <- seq(-1, 1, length.out = 500)

# How many frequencies to show (k = 1,...,K)
K <- 1  # this will give: 1, cos/sin(2πx), cos/sin(4πx), cos/sin(6πx)

# Build Fourier basis
# Basis functions: 
#   f0(x) = 1
#   f_{2k-1}(x) = cos(2π k x)
#   f_{2k}(x)   = sin(2π k x)
fourier_df <- data.frame(x = x)

# Constant term
fourier_df$f0 <- 1

# Add cos and sin terms
for (k in 1:K) {
  fourier_df[[paste0("cos", k)]] <- cos(2 * pi * k * x)
  fourier_df[[paste0("sin", k)]] <- sin(2 * pi * k * x)
}

# Convert to long format for ggplot
fourier_long <- fourier_df |>
  pivot_longer(
    cols = -x,
    names_to = "basis",
    values_to = "value"
  )

# Nice labels
fourier_long$basis <- factor(
  fourier_long$basis,
  levels = c("f0", paste0("cos", 1:K), paste0("sin", 1:K)),
  labels = c(
    "1",
    paste0("cos(2π·", 1:K, "x)"),
    paste0("sin(2π·", 1:K, "x)")
  )
)

# Plot the basis functions
p_fourier <- ggplot(fourier_long, aes(x = x, y = value, color = basis)) +
    geom_line(size = 1) +
    labs(
        x = "Input",
        y = "Basis function value",
    ) +
    theme_minimal() +
    theme(legend.position = 'none')

library(patchwork)
p_legen + p_bspline + p_fourier

fourier_long <- fourier_long %>%
    mutate(basis = as.numeric(as.factor(basis)), type = "Fourier basis")
df_legen <- df_legen %>%
    mutate(basis_order = as.factor(basis_order), type = "Legendre basis") %>%
    rename(basis = basis_order, x=variable) %>%
    dplyr::select(x, basis, value, type)

 p_basis <- rbind(fourier_long, df_legen) %>%
     ggplot(aes(x = x, y = value, color = basis)) +
     geom_line(size = 0.8) +
     scale_color_brewer(palette = "Dark2") +
     labs(x = "Input", y = "Value") +
     facet_wrap(~type) +
     #  theme_minimal() +
     theme(
         axis.title.x = element_text(face = "bold", size = 12),
         axis.title.y = element_text(face = "bold", size = 12),
         legend.position = "none",
         legend.title = element_text(size = 11),
         legend.text = element_text(size = 9),
         strip.text = element_text(face = "bold", size = 12),
         plot.title = element_text(face = "bold"),
         panel.grid.minor = element_blank()
     )
 ggsave('~/Downloads/p_basis.png', p_basis, width = 8, height = 4, dpi = 300)

 ########### sieve OLS ##########
 y <- read.csv("/Users/siyuansha/Documents/GitHub/sta496-SIMle/method/test_data.csv") %>% pull(x)
 x <- c(NA, y[1:499])
rg_df <- data.frame(y = y, x = x, t=c(1:500)) %>% data.table()
esti_df = esti_beta(y, c=5, d=5, b_time='Legen', b_timese='Legen', mp_type='algeb', r=1, s=1)
beta <- esti_df$fit_beta
W <- esti_df$d_m
est <- fix_t_esti(y, c=5, d=5, t=0, b_time='Legen', b_timese='Legen', mp_type = 'algeb', r=1, upper=1)
esti <- est$esti[[1]]
fitted_basis <- matrix(ncol = 25, nrow=2000)
for (i in 1:25) {
    basis_col <- est$basis_matrix[, i]
    fitted_basis[, i] <- beta[i] * basis_col
}
plot(est$basis_matrix[,7])

p_ols <- data.table(y = y[1:499], x = W[, 1]) %>% ggplot(aes(x = x, y = y)) +
    geom_point(color='#999999') +
    geom_smooth(method = "lm", se=FALSE) +
    labs(x = "Basis function value", y = "Y") +
    scale_color_brewer(palette = "Dark2") +
    theme(
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "none",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    ); p_ols # caption: OLS estimates coefficient for sieve basis

df <- data.table(matrix(fitted_basis, nrow = 2000, ncol = 25)) %>%
    dplyr::select(1:3) %>%
    mutate(x = seq(-1, 1, length.out = 2000), esti = V1 + V2 + V3) %>%
    melt(id.vars='x')
df

library(ggplot2)
library(dplyr)

library(ggplot2)

p_sieve_est <- ggplot(df, aes(x = x, y = value)) +
  # 1) Basis functions: V1, V2, V3 in Dark2-like colors, semi-transparent
  geom_line(
    data = subset(df, variable != "esti"),
    aes(color = variable),
    alpha = 0.4,
    linewidth = 0.7
  ) +
  # 2) Estimated function: bold black line, mapped to color = "esti"
  geom_line(
    data = subset(df, variable == "esti"),
    aes(color = "esti"),
    linewidth = 1
  ) +
  
  # Manual colors + legend labels
  scale_color_manual(
    name = "Component at t=0",  # legend title
    values = c(
      "V1"  = "#1B9E77",  # Dark2[1]
      "V2"  = "#D95F02",  # Dark2[2]
      "V3"  = "#7570B3",  # Dark2[3]
      "esti" = "black"    # estimated function
    ),
    labels = c(
      "V1"  = expression(hat(beta)[1]*b[1](t,x)),
      "V2"  = expression(hat(beta)[2]*b[2](t,x)),
      "V3"  = expression(hat(beta)[3]*b[3](t,x)),
      "esti" = 'Estimated function'
    )
  ) +
  
  labs(
    x = "Input (x)",
    y = "Value"
  ) +
  theme(
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face='bold',size = 11),
    legend.text  = element_text(size = 9),
    strip.text   = element_text(face = "bold", size = 12),
    plot.title   = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ); p_sieve_est
ggsave("~/Downloads/p_ols.png", p_ols, width = 5, height = 3)

ggsave("~/Downloads/p_sieve_est.png", p_sieve_est, width = 10, height = 3)
p_legen <- df_legen %>%  mutate(basis_num = as.integer(basis)) %>% filter(basis_num <= 3) %>%
    ggplot(aes(x = x, y = value, color = basis, alpha=0.3)) +
    geom_line(size = 0.8) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Input", y = "Value") +
    facet_wrap(~type) +
    #  theme_minimal() +
    theme(
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "none",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    ); p_legen
ggsave('~/Downloads/p_lengendre_basis.png', p_legen, width=3, height=2)

# pointwise vs. uniform
# install.packages("ggplot2") # if you don't have it
library(ggplot2)

# Suppose you have:
estimate   <- 0.8          # point estimate
se         <- 0.05         # standard error
truth      <- 0.75         # ground truth value
alpha      <- 0.05

# 95% CI
z_value <- qnorm(1 - alpha/2)  # 1.96
lower   <- estimate - z_value * se
upper   <- estimate + z_value * se

df <- data.frame(
  Method    = "x=x0",
  Estimate  = estimate,
  CI_lower  = lower,
  CI_upper  = upper,
  Truth     = truth
)


p_ci_point <-  ggplot(df, aes(x = Method, y = Estimate)) +
    
    # 95% CI error bars
    geom_errorbar(
        aes(ymin = CI_lower, ymax = CI_upper, color = "95% CI"),
        width = 0.15, linewidth = 0.7
    ) +

    # Estimation point
    geom_point(
        aes(color = "Estimation"),
        size = 3
    ) +

    # Truth point (thicker + distinct color)
    geom_point(
        aes(y = Truth, color = "Truth"),
        size = 3.5
    ) +

    # Labels
    labs(
        x = "",
        y = "f(x0)",
        color = ""
    ) +

    # Matching color scheme
    scale_color_manual(
        values = c(
            "Estimation" = "#1F1F1F",
            "95% CI"     = "#1E88E5",
            "Truth"      = "#D81B60"
        )
    ) +

    # Theme styling
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(face='bold', size=12),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "top",
        legend.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank()
    )
ggsave('~/Downloads/ci_point.png', p_ci_point, width=4, height=3)
