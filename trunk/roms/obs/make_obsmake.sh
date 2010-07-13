#!/bin/sh
set -ex
PGM=obsmake.s01
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
$F90 $OMP $F90OPT $INLINE $INCLUDE_NETCDF -c common_roms.f90
$F90 $OMP $F90OPT $INLINE -c common_obs_roms.f90
$F90 $OMP $F90OPT -c obsmake.f90
$F90 $OMP $F90OPT -o ${PGM} *.o $F90LIB $LIB_NETCDF

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

#mv *.log *.L Log/
