#!/bin/sh
set -ex
MEM=063.4d
OPT=vopt
PGM=letkf$MEM.$OPT.m01
F90=sxmpif90
OMP=
#OMP='-P openmp'
DEBUG='-C debug'
F90OPT='-C vopt -Wf"-pvctl matmulblas" -ftrace -Wf"-pvctl fullmsg -L fmtlist transform"'
#F90OPT='-C vopt -ftrace -Wf"-pvctl fullmsg -L fmtlist transform"'
#F90OPT='-C vsafe -ftrace -Wf"-pvctl fullmsg -L fmtlist transform"'
INLINE="-pi expin=common.f90"

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o
rm -rf Log$MEM.$OPT
mkdir Log$MEM.$OPT

$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT -c common_mpi.f90
$F90 $OMP $F90OPT $INLINE -c common_obs.f90
$F90 $OMP $F90OPT $INLINE -c common_mtx.f90
$F90 $OMP $F90OPT $INLINE -c netlib.f
$F90 $OMP $F90OPT -c common_letkf.f90
$F90 $OMP $F90OPT $INLINE -c common_afes.f90
$F90 $OMP $F90OPT -c common_mpi_afes.f90
$F90 $OMP $F90OPT -c common_obs_afes.f90
$F90 $OMP $F90OPT $INLINE -c letkf_tools.f90
$F90 $OMP $F90OPT -c letkf.f90
$F90 $OMP $F90OPT -o ${PGM} *.o -lblas

mv *.L Log$MEM.$OPT
rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
