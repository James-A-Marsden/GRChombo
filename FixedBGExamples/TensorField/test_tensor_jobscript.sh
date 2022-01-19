#!/bin/bash -l

#SBATCH --nodes 1

# NB cosma7 has 28 cores per node so product of #s below = 28
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=2

#SBATCH -J tensor #Give it something meaningful.

#SBATCH -o standard_output_file.%J.out
#SBATCH -e standard_error_file.%J.err

#SBATCH -p cosma7-pauper #or some other partition, e.g. cosma, cosma6, etc.

#SBATCH -A dp092 #e.g. dp004
#SBATCH -t 10:00

#Set working directory
#SBATCH -D ./

#SBATCH --mail-type=END # ALL/BEGIN/END/FAIL/NONE
#SBATCH --mail-user=james.marsden@physics.ox.ac.uk #PLEASE PUT YOUR EMAIL ADDRESS HERE (without the <>)


#load the modules used to build your program.
module purge 
module load intel_comp/2019 intel_mpi/2019 parallel_hdf5/1.10.3


# Run the program
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np $SLURM_NTASKS ./Main_TensorField3d.Linux.64.mpiicpc.ifort.OPT.MPI.OPENMPCC.ex params_tensor_test.txt
