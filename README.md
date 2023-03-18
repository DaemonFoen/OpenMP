# OpenMP
A program for solving a system of linear algebraic equations using OpenMP

qsub script:

#!/bin/sh
#PBS -l walltime=00:00:50
#PBS -l select=1:ncpus=12:ompthreads=1
#PBS -e err.txt
#PBS -o log.txt

cd $PBS_O_WORKDIR
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo
./omp_one
