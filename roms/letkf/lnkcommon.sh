#!/bin/sh
# for making link of common sources
set -e
COMMONDIR=../../common
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/common_mtx.f90 ./
ln -fs $COMMONDIR/common_letkf.f90 ./
ln -fs $COMMONDIR/netlib.f ./

ln -fs ../common/common_roms.f90 ./
ln -fs ../common/common_mpi_roms.f90 ./
ln -fs ../common/common_obs_roms.f90 ./

