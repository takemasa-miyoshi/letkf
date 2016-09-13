#!/bin/bash
#PBS -l nodes=4:ppn=2

export OMP_NUM_THREADS=1
export PARALLEL=1

export FORT_FMT_RECL=400

. /usr/share/modules/init/sh
module unload pgi-12.10
module load intel-2013.1.117
module load mpt/2.10
module load hdf5/1.8.13
module load netcdf/4.3.2
module load netcdf-fortran/4.2

ulimit -s unlimited
umask 022

cd $PBS_O_WORKDIR

rm -f machinefile
cp -f $PBS_NODEFILE machinefile

cat machinefile | sort | uniq > mpd.hosts
nnodes=`cat mpd.hosts | wc -l`

mpdboot -n $nnodes -f mpd.hosts
#mpdtrace

nodelist=''
for inode in `cat mpd.hosts`; do
  if [ -z "$nodelist" ]; then
    nodelist="$inode"
  else
    nodelist="$nodelist,$inode"
  fi
done

export MPI_UNIVERSE="$nodelist 40"
echo "MPI_UNIVERSE='$MPI_UNIVERSE'"

#-------------------------------------

./fcst.sh

#-------------------------------------

mpdallexit

exit 0
