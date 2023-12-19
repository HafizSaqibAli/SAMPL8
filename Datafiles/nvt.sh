#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N nvt
#$ -pe smp.pe 8  
 
module load apps/intel-17.0/gromacs/2018.4/single

export OMP_NUM_THREADS=$NSLOTS
gmx mdrun -deffnm nvt

