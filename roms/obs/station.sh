#!/bin/sh
set -ex
PGM=station.s01
F90=ifort
OMP=
F90OPT='-O3'
INLINE=
#INCLUDE_NETCDF="-I/usr/local/netcdf/include"
#LIB_NETCDF="-L/usr/local/netcdf/lib -lnetcdf"
INCLUDE_NETCDF="-I/home/zhc/netcdf-3.6.0-p1/include"
LIB_NETCDF="-L/home/zhc/netcdf-3.6.0-p1/lib -lnetcdf"

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
