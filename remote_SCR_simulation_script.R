path_to_helper <- '/scratch/kssha/sta496/SCR_sim_helper.R'
save_path <- '/scratch/kssha/sta496/sim_res/'
source(path_to_helper)

args <- commandArgs(trailingOnly = TRUE)
n_rep <- as.integer(args[1])
n <- as.integer(args[2])
m_case <- as.integer(args[3])
eps_case <- args[4]
seed <- as.integer(args[5])
n_cores <- as.integer(args[6])  # New: get cores from SLURM
cat('Running with n_rep=', n_rep, ' n=', n, ' m_case=', m_case, ' eps_case=', eps_case, ' seed=', seed, ' n_cores=', n_cores, '\n')
c=7
d=7

n_options <- 50 # c(500, 800)
eps_options <- c('tvAR2', 'SETAR', 'bilinear')

cov <- auto_SCR_coverage(n_rep = n_rep, n=n, c=c, d=d, delta=1, m_case = m_case, eps_case = eps_case, n_cores = n_cores, debug=TRUE)
file_name <- paste0(save_path, paste0('n_reps=', n_rep,'n=',n,'c=',c,'d=',d, 'm', m_case, 'eps=', eps_case,
Sys.Date(), '.rds'))

saveRDS(cov, file_name)
cat('Saved to ', file_name)

# for (n in n_options) {
#   res <- list()
#   
#   for (eps in eps_options) {
#     cat('Running case n=', n, eps, '\n')
#     cov <- auto_SCR_coverage(n_rep = n_rep, n=n, f_true=m2, delta=1, m_case = m_case, eps_case = eps)
#     res[paste0('n=', n,'_', eps)] <- cov
#   }
#   
#   file_name <- paste0(save_path, paste0('n=', n, 'n_reps=', n_rep, 'm', m_case, 'eps=', paste(eps_options, collapse = '-'), '.rds'))
#   
#   saveRDS(res, file_name)
#   cat('Saved to ', file_name)
# }
