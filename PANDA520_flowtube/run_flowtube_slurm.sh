#!/bin/bash
# Example SLURM job script for running MARFORCE-Flowtube on an HPC cluster.
# Edit the account, module and path lines marked "EDIT" for your system.
#
# Usage (from PANDA520_flowtube/):
#   mkdir -p logs
#   sbatch run_flowtube_slurm.sh Start_SetParam_Isoprene.py
#
#SBATCH --account=project_XXXXXXX      # EDIT: your project/account
#SBATCH --job-name=flowtube
#SBATCH --partition=small
#SBATCH --time=03-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --output=logs/log_%x_%j.out

# ── Environment (EDIT) ──
# Either load your cluster's Python + Open Babel modules and a virtualenv ...
#   module load python-data/3.10-24.04
#   module load openbabel
#   source /path/to/my_env/bin/activate
# ... or activate the conda environment created from environment.yml:
#   source ~/miniconda3/etc/profile.d/conda.sh && conda activate marforce

# ── Thread control ──
# The flow-tube solver parallelises over grid cells with Numba; keep BLAS single-threaded.
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR:-$PWD/.numba_cache}
mkdir -p "$NUMBA_CACHE_DIR"
export MPLBACKEND=Agg   # no display on compute nodes

# ── Run ──
cd "$SLURM_SUBMIT_DIR"
export PYTHONPATH="$(pwd):$PYTHONPATH"
SCRIPT=${1:-Start_SetParam_Isoprene.py}
echo "Running: $SCRIPT"
echo "Job ID: $SLURM_JOB_ID   CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
python "$SCRIPT"
echo "End: $(date)"
