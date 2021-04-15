#!/bin/bash

#SBATCH --nodes=6
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --output=mpi-out.%j
#SBATCH --error=mpi-err.%j
#SBATCH --time=01:30:00
#SBATCH --partition=pool

rm -f gradients/grad_checkerboard_dh10cm_Surface_Stream_Elastic2D*
rm -f gradients/Hessian_checkerboard_dh10cm_Surface_Stream_Elastic2D*
rm -f seismograms/seismograms_checkerboard_dh10cm_Surface_Stream_Elastic2D_Start.stage_*
rm -f seismograms/seismograms_checkerboard_dh10cm_Surface_Stream_Elastic2D_True.stage_*
rm -f model/model_checkerboard_dh10cm_Surface_Stream_Elastic2D_Start.stage_*
rm -f sourceInversion/invSource_checkerboard_dh10cm_Surface_Stream_Elastic2D*
rm -f logs/steplengthSearch_checkerboard_dh10cm_Surface_Stream_Elastic2D.log

export OMP_NUM_THREADS=1
export SCAI_UNSUPPORTED=IGNORE
export SCAI_TRACE=OFF

mpirun -np 1 ./../build/bin/Inversion "configuration/configuration_checkerboard_dh10cm_Surface_Stream_Elastic2D_Inv_HPC.txt"
