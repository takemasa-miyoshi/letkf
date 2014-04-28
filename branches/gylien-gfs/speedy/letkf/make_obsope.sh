#!/bin/sh
set -ex
PGM=obsope.s01
F90=pgf90
OMP=
F90OPT='-byteswapio -tp sandybridge-64 -fast -O3'
INLINE="-Minline"
BLAS=1 #0: no blas 1: using blas

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT $INLINE -c common_speedy.f90
$F90 $OMP $F90OPT -c common_obs_speedy.f90
$F90 $OMP $F90OPT -c obsope.f90
$F90 $OMP $F90OPT -o ${PGM} *.o

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
