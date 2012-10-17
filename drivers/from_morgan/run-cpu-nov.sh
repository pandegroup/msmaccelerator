#!/bin/bash

#PBS -N gyrB_nov_run10
#PBS -e run10-nov.err
#PBS -o run10-nov.out
#PBS -l nodes=1:ppn=24
#PBS -M aparente@stanford.edu 
#PBS -m e 
#PBS -l walltime=24:00:00
#PBS -V 

#export MYNUMPROCS=8
export AMBERHOME=/home/mlawrenz/Programs/amber12/bin/
echo "USING $MYNUMPROCS processors for runs\n"
export OMP_NUM_THREADS=12

WORKDIR=$WORKDIR/

#Minimization
/share/apps/openmpi-1.4.2/intel-11/bin/mpirun -machinefile  $PBS_NODEFILE -np 24 /home/mlawrenz/Programs/amber12/bin/pmemd.MPI -O -i $WORKDIR/min-nov.in -o $WORKDIR/min-nov.out -p $WORKDIR/aln-4dua-model-complex.solv.prmtop -c $WORKDIR/aln-4dua-model-complex.solv.inpcrd -r $WORKDIR/min-nov.rst -ref $WORKDIR/aln-4dua-model-complex.solv.inpcrd 
#Heating
/share/apps/openmpi-1.4.2/intel-11/bin/mpirun -machinefile $PBS_NODEFILE -np 24 /home/mlawrenz/Programs/amber12/bin/pmemd.MPI  -O -i $WORKDIR/heat-nov.in -o $WORKDIR/heat-nov.out -p $WORKDIR/aln-4dua-model-complex.solv.prmtop -c $WORKDIR/min-nov.rst -r $WORKDIR/heat-nov.rst -x $WORKDIR/heat-nov.mdcrd -ref $WORKDIR/min-nov.rst  
#Density
/share/apps/openmpi-1.4.2/intel-11/bin/mpirun -machinefile $PBS_NODEFILE -np 24 /home/mlawrenz/Programs/amber12/bin/pmemd.MPI -O -i $WORKDIR/density-nov.in -o $WORKDIR/density-nov.out -p $WORKDIR/aln-4dua-model-complex.solv.prmtop -c $WORKDIR/heat-nov.rst -r $WORKDIR/density-nov.rst -x $WORKDIR/density-nov.mdcrd -ref $WORKDIR/heat-nov.rst
#Equilibration
/share/apps/openmpi-1.4.2/intel-11/bin/mpirun -machinefile $PBS_NODEFILE -np 24 /home/mlawrenz/Programs/amber12/bin/pmemd.MPI -O -i $WORKDIR/equil-nov.in -o $WORKDIR/equil-nov.out -p $WORKDIR/aln-4dua-model-complex.solv.prmtop -c $WORKDIR/density-nov.rst -r $WORKDIR/equil-nov.rst -x $WORKDIR/equil-nov.mdcrd

