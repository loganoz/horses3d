#!/bin/bash
#SBATCH --job-name=iLES
#SBATCH --nodelist=n002
#SBATCH --ntasks=20
#SBATCH --time=24:00:00
#SBATCH -e RESULTS/err%j.log
#SBATCH -o RESULTS/out%j.log
 
export OMP_SCHEDULE="guided"
export OMP_NUM_THREADS=1
 
module purge
module load gcc/11.2.0
module load hdf5/gcc11.2/1.14.2
module load openmpi/gcc11.2/4.1.6
module load metis/gcc11.2/5.1.0
 
echo "##########################################################################"
echo "#"
echo "# Running with $SLURM_NTASKS tasks and $SLURM_CPUS_PER_TASK cpus/task"
echo "# On nodes $SLURM_JOB_NODELIST"
echo "#"
echo "##########################################################################"
 
EXEC=../../../bin/horses3d.ns
 
mpiexec $EXEC Cylinder.control
