#!/bin/sh
set -ex
INPUT=ee6_rst.0000.nc
F90=pgf90
PGM=nc2grd
LIB_NETCDF="-L/usr/local/netcdf/lib -lnetcdf"
INC_NETCDF="-I/usr/local/netcdf/include"
ln -fs ../../common/SFMT.f90 .
ln -fs ../../common/common.f90 .
ln -fs ../common/common_roms.f90 .
F90OPT='-byteswapio -tp p7-64 -fast -O3'
$F90 $F90OPT -c SFMT.f90
$F90 $F90OPT -c common.f90
$F90 $F90OPT $INC_NETCDF -c common_roms.f90
$F90 $F90OPT $INC_NETCDF -c nc2grd.f90
$F90 -o $PGM *.o $LIB_NETCDF
rm SFMT.f90
rm common.f90
rm common_roms.f90
rm *.mod
rm *.o
ln -fs $INPUT input.nc
cp $INPUT output.nc
./$PGM
rm input.nc
rm $PGM
