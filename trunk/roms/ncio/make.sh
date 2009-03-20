#!/bin/sh
set -ex
F90=pgf90
#PGM=chtimestep
PGM=choceantime
LIB_NETCDF="-L/usr/local/netcdf/lib -lnetcdf"
INC_NETCDF="-I/usr/local/netcdf/include"
F90OPT='-byteswapio -tp p7-64 -fast -O3'
$F90 $F90OPT $INC_NETCDF -o $PGM $PGM.f90 $LIB_NETCDF
rm *.o
