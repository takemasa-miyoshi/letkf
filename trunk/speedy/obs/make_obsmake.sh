#!/bin/sh
set -ex
PGM=obsmake.s01
F90=pgf90
OMP=
F90OPT='-byteswapio -tp p7-64 -fast -O3'
INLINE="-Minline"

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT $INLINE -c common_speedy.f90
$F90 $OMP $F90OPT -c common_obs_speedy.f90
$F90 $OMP $F90OPT -c obsmake.f90
$F90 $OMP $F90OPT -o ${PGM} *.o $F90LIB

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
