#!/bin/sh
F90=ifort
OBS=../DATA/obs/regular/2004010100_obs.dat
ln -s $OBS fort.3
$F90 obsdump.f90
./a.out
rm a.out
rm fort.3

