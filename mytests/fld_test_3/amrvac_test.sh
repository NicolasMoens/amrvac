#!/bin/sh
#SBATCH --job-name=amrvac_test
#SBATCH --account=ivsusers
#SBATCH --time 300
#  (estimated run time in minutes)
#SBATCH --tasks-per-node=1
#SBATCH -N 1
#  (default value, use this value if you want to execute a job on multiple nodes (openmpi))
#SBATCH --mem=2048
#  (memory in MB)
#SBATCH --partition=normal
#   (use the normal partition=default)
#SBATCH --output=/home/nicolasm/amrvac/mytests/fld_test_3/stdout.log
#SBATCH --error=/home/nicolasm/amrvac/mytests/fld_test_3/stderr.log


export AMRVAC_DIR=$HOME/amrvac 
export PATH="/home/nicolasm/MPICH_3.2/bin:$PATH"

cd /home/nicolasm/amrvac/mytests/fld_test_3 
mpirun /home/nicolasm/amrvac/mytests/fld_test_3/amrvac -i /home/nicolasm/amrvac/mytests/fld_test_3/usr2.par
srun date
