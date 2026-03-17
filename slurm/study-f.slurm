#!/bin/bash
#SBATCH -J reg_study_f
#SBATCH -o ./logs/%x.%j.out
#SBATCH -e ./logs/%x.%j.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/sim-registration
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=48000mb
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE

module load slurm_setup
module load r/4.3.3-gcc13-mkl

export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.3"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

mkdir -p logs results
Rscript sim-run.R F 17
