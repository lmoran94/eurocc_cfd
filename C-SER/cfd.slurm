#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=cfd
#SBATCH --time=0:50:0
#SBATCH --ntasks=1
#SBATCH --account=dc139 
#SBATCH --partition=standard
#SBATCH --qos=standard

# Load the default HPE MPI environment
module load intel-compilers-19
module load mpt

# Change to the submission directory
cd $SLURM_SUBMIT_DIR

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.

for t in 1 2 4 8 16 32
do
	export OMP_NUM_THREADS=$t
 
# Launch the serial job
#   Using 1 thread
	srun --cpu-bind=cores ./cfd 32 5000
done
