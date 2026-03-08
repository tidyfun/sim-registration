#!/bin/bash
#SBATCH -J tf_fda_rerun
#SBATCH -o ./logs/%x.%j.out
#SBATCH -e ./logs/%x.%j.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/tf/attic/sim-registration
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16000mb
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE

module load slurm_setup
module load r/4.3.3-gcc13-mkl

export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

Rscript study-a-fda-rerun.R 8
