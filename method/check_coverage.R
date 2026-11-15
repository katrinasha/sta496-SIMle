library(Sie2nts)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(parallel)
library(snow)
# source('/Users/siyuansha/Documents/GitHub/sta496-SIMle/R/SIMle.Legen.v1.R')
# source("/Users/siyuansha/Documents/GitHub/sta496-SIMle/R/SIMle.original_code.v1.R")
# source('R/fix.SCR.R')
# source('/Users/siyuansha/Documents/GitHub/sta496-SIMle/R/auto.SCR.R')
source("/Users/siyuansha/Documents/GitHub/sta496-SIMle/method/debug_helper.R")
# source('/Users/siyuansha/Documents/GitHub/sta496-SIMle/method/SIMle.original_copy.R')
library(SIMle)

m_case = 1
eps_case='tvAR2'
n = 500
c = 5
d = 5
b_time = "Legen"
b_timese = "Legen"
s = 1
r = 1
n_point = 2000
upper = 0.5
n_rep = 10
ops = "CV"
per = 0.1
mp_type='algeb'
library(SIMle)
cov_p <- auto_SCR_coverage(n_rep, n, c, d,
    delta = NA, b_time, b_timese, mp_type,
    type = "fixt", m_case = m_case, eps_case = eps_case
)
cov_p