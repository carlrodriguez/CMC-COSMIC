#!/bin/bash

#SBATCH -J Omega_Cen
#SBATCH --error=error.out
#SBATCH --output=output.out
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=300G
#SBATCH --time 700:00:00 
#SBATCH --partition=carlrlab


module purge all
module load cmake/3.18.0 
module load openmpi/4.1.4-intel_20.2

##mpirun -np X(=n_nodes*n_cores) <exe> > output 
ulimit -s unlimited
mpirun -np 128 ./cmc ElsonProfile_OmegaCen_bhlc.ini initial

