library(rvest)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(SIMle)
setwd('~/Desktop/UofT/STA496 Time Series/sta496-coding')

data <- read.csv('hkepd_dat.csv') %>% data.table()

# model setup
y <- data %>% pull(Total) %>% scale(center = TRUE, scale = FALSE)
x1 <- data %>% pull(SO2) %>% scale()
x2 <- data %>% pull(NO2) %>% scale()
x3 <- data %>% pull(RSP) %>% scale()

# time homogenous test
res <- homo.test(x1, 5, 2, m = 'MV', 'Legen', 'Legen', 'algeb', fix_num = 0, r = 2, n_point = 1000)
fix_nums <- seq(from=-10, to=10, length.out = 5)
res_df <- data.frame()
for (fix_num in fix_nums) {
  cat('fixing time for: ', fix_num)
  res <- homo.test(x1, 5, 2, m = 'MV', 'Legen', 'Legen', 'algeb', fix_num = fix_num, r = 1, n_point = 1000)[[1]] %>%
    mutate(fix_t = fix_num)
  res_df <- rbind(res_df, res)
}
head(res_df)

library(plotly)

ts_df %>%
  plot_ly(
    x = ~x_tm1,
    y = ~t,
    z = ~xt,
    type = "surface",
    mode = "markers",
    marker = list(size = 3, color = ~y, colorscale = "Viridis")
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "x_tm1"),
      yaxis = list(title = "t"),
      zaxis = list(title = "xt")
    )
  )
library(plotly)

plot_ly() %>%
  add_trace(
    data = subset(res_df, order == "esti"),
    x = ~x, y = ~fix_t, z = ~y,
    type = "scatter3d", mode = "markers",
    marker = list(size = 3, color = ~y, colorscale = "Viridis"),
    name = "Estimate"
  ) %>%
  add_trace(
    data = subset(res_df, order == "upper"),
    x = ~x, y = ~fix_t, z = ~y,
    type = "scatter3d", mode = "markers",
    marker = list(size = 3, color = ~y, colorscale = "Blues"),
    name = "Upper CI"
  ) %>%
  add_trace(
    data = subset(res_df, order == "lower"),
    x = ~x, y = ~fix_t, z = ~y,
    type = "scatter3d", mode = "markers",
    marker = list(size = 3, color = ~y, colorscale = "Reds"),
    name = "Lower CI"
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Fix_t"),
      zaxis = list(title = "Y")
    ),
    legend = list(title = list(text = "Type"))
  )
test.plot(res[[1]], 'homot')
test.plot(res[[2]], 'homot')
ts  = generate_nAR1(4000, 1)
getAnywhere(homo.test)
SIMle:::SCR_fix_t
?homo.test

