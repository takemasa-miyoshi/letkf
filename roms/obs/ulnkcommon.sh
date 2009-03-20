#!/bin/sh
# for deleting links of common sources
set -e

rm -f SFMT.f90
rm -f common.f90
rm -f common_mpi.f90
rm -f common_obs.f90
rm -f common_mtx.f90
rm -f common_letkf.f90
rm -f netlib.f

rm -f common_roms.f90
rm -f common_mpi_roms.f90
rm -f common_obs_roms.f90

