#!/bin/sh
set -ex
#PGM=letkf020.m01
PGM=test
F90=ifort
MPI='-lmpi'
OMP=
# F90OPT='-byteswapio -tp p7-64 -fast -O3'
F90OPT='-O3'
# F90OPT='-64 -parallel -omp -pardiag=3 -O4 -pmfunc -pmpar -fullmsg -loglist'
# F90OPT='-64 -parallel -omp -pardiag=3 -O4 -loglist -debug'
# F90OPT='-64 -parallel -omp -pardiag=3 -O4 -fullmsg -loglist -oplist'
# F90OPT='-64 -O4 -loopdistribute -loopexpand -loopfuse -loopinterchange -loopreroll -parallel=1 -omp -pardiag=3 -fullmsg -loglist -oplist'
#F90OPT='-64 -O4 -loopdiag -loopdistribute -loopinterchange -parallel=1 -omp -pardiag=3 -fullmsg -loglist -oplist'

INLINE=
#F90LIB="-L/users/npd/suuchi54/dvlK2/rttov87/src -lrttov8.7"
#F90LIB="-L/users/npd/suuchi54/dvlK2/rttov87_new/src -lrttov8.7"
#INCLUDE="-I/users/npd/suuchi54/dvlK2/rttov87_new/src"
#INCLUDE_NETCDF="-I/usr/local/netcdf/include"
#LIB_NETCDF="-L/usr/local/netcdf/lib -lnetcdf"
INCLUDE_NETCDF="-I/home/zhc/netcdf-3.6.0-p1/include"
LIB_NETCDF="-L/home/zhc/netcdf-3.6.0-p1/lib -lnetcdf"

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INLINE -c SFMT.f90 $MPI
$F90 $OMP $F90OPT $INLINE -c common.f90 $MPI
$F90 $OMP $F90OPT -c common_mpi.f90 $MPI
$F90 $OMP $F90OPT $INLINE -c common_mtx.f90 $MPI
$F90 $OMP $F90OPT $INLINE -c netlib.f $MPI
$F90 $OMP $F90OPT -c common_letkf.f90 $MPI
$F90 $OMP $F90OPT $INLINE $INCLUDE_NETCDF -c common_roms.f90 $MPI
$F90 $OMP $F90OPT -c common_obs_roms.f90 $MPI
$F90 $OMP $F90OPT -c common_mpi_roms.f90 $MPI
$F90 $OMP $F90OPT -c letkf_obs.f90 $MPI
$F90 $OMP $F90OPT -c letkf_tools.f90 $MPI
$F90 $OMP $F90OPT -c letkf.f90 $MPI
$F90 $OMP $F90OPT -o ${PGM} *.o $F90LIB $LIB_NETCDF $MPI

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

#mv *.log *.L Log/
