#!/bin/sh
set -ex
PGM=station.s01
F90=pgf90
OMP=
F90OPT='-byteswapio -tp p7-64 -fast -O3'
INLINE="-Minline"
INCLUDE_NETCDF="-I/usr/local/netcdf/include"
LIB_NETCDF="-L/usr/local/netcdf/lib -lnetcdf"

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT $INLINE -c common_obs.f90
$F90 $OMP $F90OPT $INLINE $INCLUDE_NETCDF -c common_roms.f90
$F90 $OMP $F90OPT $INCLUDE -c station.f90
$F90 $OMP $F90OPT -o ${PGM} *.o $F90LIB $LIB_NETCDF

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh
### run
./$PGM > station.tbl
rm $PGM
