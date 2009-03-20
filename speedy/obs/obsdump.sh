#!/bin/sh
F90=pgf90
OBS=../DATA/obs/1982010100.dat
ln -s $OBS fort.3
$F90 -byteswapio obsdump.f90
./a.out
rm a.out
rm fort.3

