#!/bin/bash
#SBATCH --account=def-zhouzhou        # replace this with your own account
#SBATCH --time=03:00:00
#SBATCH --job-name="SCR simulation" # set to your job name
#SBATCH --nodes=1                   # keep to 1 for shared-memory parallel
#SBATCH --cpus-per-task=192         # number of CPU cores per task (adjust based on your needs)
#SBATCH --mail-user=katrina.sha@mail.utoronto.ca
#SBATCH --mail-type=ALL 
#SBATCH --output=debug/%x_%A_%a.out
#SBATCH --error=debug/%x_%A_%a.err
#SBATCH --array=1-9                # product of grids (1*1*3*3)

# Increase maximum number of open file descriptors
ulimit -n 65536

# Set OpenMP/MKL/BLAS thread count to 1 (important for R parallel)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

module load StdEnv/2023 r/4.5.0
cd /scratch/kssha/sta496
export R_LIBS=~/local/R_libs/

# Define your grids
N_REP_LIST=(300)
N_LIST=(100)

# N_LIST can be provided externally via the environment variable N_LIST_ARG.
# Accepts either a comma-separated string (e.g. "100,500") or a space-separated
# string (e.g. "100 500"). If not provided, defaults to (500).
if [ -n "${N_LIST_ARG:-}" ]; then
  if [[ "$N_LIST_ARG" == *,* ]]; then
    IFS=',' read -r -a N_LIST <<< "$N_LIST_ARG"
  else
    read -r -a N_LIST <<< "$N_LIST_ARG"
  fi
  echo "N_LIST overridden from environment: ${N_LIST[*]}"
else
  N_LIST=(500)
  echo "Using default N_LIST: ${N_LIST[*]}"
fi

M_LIST=(1 2 3)
EPS_LIST=(tvAR2 SETAR bilinear)

# Compute indices from SLURM_ARRAY_TASK_ID (0-based math)
IDX=$((SLURM_ARRAY_TASK_ID - 1))

a=$(( IDX % ${#N_REP_LIST[@]} ))
IDX=$(( IDX / ${#N_REP_LIST[@]} ))

b=$(( IDX % ${#N_LIST[@]} ))
IDX=$(( IDX / ${#N_LIST[@]} ))

c=$(( IDX % ${#M_LIST[@]} ))
IDX=$(( IDX / ${#M_LIST[@]} ))

d=$(( IDX % ${#EPS_LIST[@]} ))

N_REP=${N_REP_LIST[$a]}
N=${N_LIST[$b]}
M_CASE=${M_LIST[$c]}
EPS_CASE=${EPS_LIST[$d]}

SEED=$((2000 + SLURM_ARRAY_TASK_ID))
echo "Task $SLURM_ARRAY_TASK_ID -> n_rep=$N_REP n=$N m_case=$M_CASE eps_case=$EPS_CASE seed=$SEED"

echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
echo "Running on host: $(hostname)"
echo "Available CPU cores: $(nproc)"

R CMD BATCH --no-save --no-restore \
  "--args $N_REP $N $M_CASE $EPS_CASE $SEED $SLURM_CPUS_PER_TASK" \
  SCR_simulation_script.R \
  debug/run_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.Rout
  
# sbatch SCR_sim_array.sh
