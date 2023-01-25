#!/bin/bash

#PBS -L tasks=4:lprocs=40:gpus=1:default
#PBS -L tasks=4:lprocs=40:gpus=4:default
#PBS -L tasks=8:lprocs=24:gpus=2:default:feature=teslav100
#PBS -q coc-ice-gpu
#PBS -j oe

cd $PBS_O_WORKDIR
cat $PBS_GPUFILE
awk -F. '{slots[$1]+=1} END {for(host in slots) print host" slots="slots[host]}' $PBS_GPUFILE > hosts

#module loads here
module load gcc
module load openmpi/4.1.2
module load nvhpc
module load cuda/11.2
module load openmpi/4.1.2

make clean
make

mpirun -np $(wc -l <${PBS_GPUFILE}) -bind-to none --hostfile hosts -x LD_LIBRARY_PATH ./run
