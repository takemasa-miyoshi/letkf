#!/bin/sh
set -e
F90=pgf90
#ORO='_oro'
ORO=
CDIR=`pwd`
cd ../..
L96DIR=`pwd`
cd ..
ENKFDIR=`pwd`
COMDIR=$ENKFDIR/common
OUTDIR=$L96DIR/DATA
WKDIR=$L96DIR/tmp
rm -rf $WKDIR
mkdir -p $WKDIR
cd $WKDIR
cp $COMDIR/SFMT.f90 .
cp $COMDIR/common.f90 .
cp $L96DIR/model/lorenz96$ORO.f90 .
cp $L96DIR/model/run/nature.f90 .
$F90 -o nature SFMT.f90 common.f90 lorenz96$ORO.f90 nature.f90
rm *.mod
rm *.o
ln -s $OUTDIR/spinup/init.dat fort.10
time ./nature
mv fort.90 $OUTDIR/nature.dat
cp $CDIR/nature.ctl $OUTDIR
