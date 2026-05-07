#!/bin/bash
#SBATCH --job-name=$SLURM_JOB_NAME
#SBATCH --partition=standard
#SBATCH --nodes=1
##SBATCH --nodelist=n006
#SBATCH  --ntasks-per-node=4
#SBATCH --mem=20G
# #SBATCH --time=24:00:00
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.capel@upm.es
#SBATCH -e RESULTS/err-al%j.log
#SBATCH -o RESULTS/out-al%j.log
 
export OMP_SCHEDULE="guided"
export OMP_NUM_THREADS=1
 
module purge
module load intel/2021a
module load HDF5/1.10.7-iimpi-2021a
module load METIS/5.1.0-foss-2018b
export HDF5_ROOT=$HDF5_DIR
export METIS_HOME=/sw/software/METIS/5.1.0-foss-2018b
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
 
echo "##########################################################################"
echo "#"
echo "# Running with $SLURM_NTASKS tasks and $SLURM_CPUS_PER_TASK cpus/task"
echo "# On nodes $SLURM_JOB_NODELIST"
echo "#"
echo "##########################################################################"
 
# EXEC=../../bin/horses3d.ns 
# ## srun $EXEC WPFS_wing_pylon_v2.control
# srun -n $SLURM_NTASKS $EXEC WPFS_moving_store_IBM.control

EXEC=../../bin/horses2plt
srun $EXEC Postpro_volume.control
srun $EXEC Postpro_surface.control