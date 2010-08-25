#!/bin/sh
set -ex
F90=pgf90
F90OPT='-byteswapio -tp p7-64 -fast -O3'
LIB_NETCDF="-L/usr/local/netcdf/lib -lnetcdf"
INC_NETCDF="-I/usr/local/netcdf/include"
###
### nc2grd
###
PGM=nc2grd
test -f $PGM && rm $PGM
$F90 $F90OPT $INC_NETCDF -c nc2grd.f90
$F90 -o $PGM *.o $LIB_NETCDF
rm *.o
###
### init_merge
###
PGM=init_merge
test -f $PGM && rm $PGM
$F90 $F90OPT $INC_NETCDF -c init_merge.f90
$F90 -o $PGM *.o $LIB_NETCDF
rm *.o
echo "NORMAL END"
