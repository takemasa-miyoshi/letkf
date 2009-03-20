#!/bin/sh
F90=pgf90
$F90 station.f90
./a.out > station.tbl
rm a.out
