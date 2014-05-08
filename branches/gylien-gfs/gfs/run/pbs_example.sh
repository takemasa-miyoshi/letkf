#!/bin/bash

#PBS -S /bin/bash
#PBS -N gfs-letkf
#PBS -l walltime=12:00:00
#PBS -l select=8:ncpus=16:mpiprocs=16
#PBS -q general
#PBS -W group_list=s1314
#PBS -W umask=027

cd $PBS_O_WORKDIR

rm -f machinefile
cp -f $PBS_NODEFILE machinefile
cat machinefile | sort | uniq > mpd.hosts
nnodes=`cat mpd.hosts | wc -l`

mpdboot -n $nnodes -f mpd.hosts
#mpdtrace


./mcycle.sh 2008010100 2008010118


mpdallexit

exit 0
