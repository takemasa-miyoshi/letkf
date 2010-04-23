#!/bin/sh
set -ex
PGM=letkf020.m01
F90=mpif90
OMP=
F90OPT='-byteswapio -tp p7-64 -fast -O3'
INLINE="-Minline"

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT -c common_mpi.f90
$F90 $OMP $F90OPT $INLINE -c common_mtx.f90
$F90 $OMP $F90OPT $INLINE -c netlib.f
$F90 $OMP $F90OPT -c common_letkf.f90
$F90 $OMP $F90OPT $INLINE -c common_speedy.f90
$F90 $OMP $F90OPT -c common_obs_speedy.f90
$F90 $OMP $F90OPT -c common_mpi_speedy.f90
$F90 $OMP $F90OPT -c letkf_obs.f90
$F90 $OMP $F90OPT -c letkf_tools.f90
$F90 $OMP $F90OPT -c letkf.f90
$F90 $OMP $F90OPT -o ${PGM} *.o

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
