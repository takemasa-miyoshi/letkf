#!/bin/bash
#PBS -l select=36:ncpus=1

export FORT_FMT_RECL=400

. /usr/share/modules/init/sh
module unload pgi-12.10
module load intel-2013.1.117
module load mpt/2.10
module load hdf5/1.8.13
module load netcdf/4.3.2
module load netcdf-fortran/4.2

ulimit -s unlimited

cd $PBS_O_WORKDIR

time mpiexec_mpt -n 36 dplace -s1 scale-rm run.conf
