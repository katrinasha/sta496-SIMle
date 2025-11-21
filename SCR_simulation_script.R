path_to_helper <- '/scratch/kssha/sta496/SCR_sim_helper.R'
save_path <- '/scratch/kssha/sta496/sim_res/'
source(path_to_helper)

cat("\nR Script Starting\n==================\n")

# Get and show all arguments
args <- commandArgs(trailingOnly = TRUE)
cat("Raw arguments received:", paste(args, collapse=" "), "\n\n")

# Show all command line arguments for debugging
all_args <- commandArgs(trailingOnly = FALSE)
cat("All R command line arguments:\n")
print(all_args)
cat("\n")

# Validate number of arguments
if (length(args) < 6) {
    stop("Not enough arguments. Expected 6 args: n_rep n m_case eps_case seed n_cores\n",
         "Received: ", length(args), " args: ", paste(args, collapse=" "))
}

# Parse arguments with validation
n_rep <- as.integer(args[1])
n <- as.integer(args[2])
m_case <- as.integer(args[3])
eps_case <- args[4]
seed <- as.integer(args[5])
n_cores <- as.integer(args[6])

# Validate critical values
stopifnot(!is.na(n_rep), !is.na(n), !is.na(m_case), !is.na(seed))

# If n_cores is NA, use a safe default
if (is.na(n_cores)) {
    warning("n_cores was NA, defaulting to min(detectCores()-1, n_rep)")
    n_cores <- min(parallel::detectCores() - 1, n_rep)
}

cat("\nValidated arguments:\n",
    "n_rep    = ", n_rep, "\n",
    "n        = ", n, "\n",
    "m_case   = ", m_case, "\n",
    "eps_case = ", eps_case, "\n",
    "seed     = ", seed, "\n",
    "n_cores  = ", n_cores, "\n",
    "Running on: ", Sys.info()["nodename"], "\n",
    "SLURM_CPUS_PER_TASK = ", Sys.getenv("SLURM_CPUS_PER_TASK"), "\n",
    sep="")

n_options <- 50 # c(500, 800)
eps_options <- c('tvAR2', 'SETAR', 'bilinear')

cov <- auto_SCR_coverage(n_rep = n_rep, n=n, delta=1, m_case = m_case, 
eps_case = eps_case, n_cores = n_cores, debug=TRUE)

auto_SCR_coverage(n_rep=1, n=500, delta=1, m_case = 1, eps_case = 'tvAR2')

file_name <- paste0(save_path, paste0('n_reps=', n_rep,'n=',n, 'm', m_case, 'eps=', eps_case, '.rds'))

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
