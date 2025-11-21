#!/bin/bash
#SBATCH --account=def-jgrons   # replace this with your own account
#SBATCH --time=00:15:00
#SBATCH --job-name="SCR simulation"  # set to your job name
#SBATCH --nodes=1  # number of node to call
#SBATCH --mail-user=katrina.sha@mail.utoronto.ca
#SBATCH --ntasks=40  # number of jobs to be run in each node
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load StdEnv/2023 r/4.5.0
cd /scratch/kssha/sta496
export R_LIBS=~/local/R_libs/

R CMD BATCH --no-save "--args $n_rep $n $m_case $eps_case" /scratch/kssha/sta496/SCR_simulation_script.R /scratch/kssha/sta496/SCR_simulation_script.Rout
# sbatch --export=n_rep=1,n=10,m_case=1,eps_case='tvAR2' SCR_sim.sh