#!/bin/sh
set -e
ORO=
OBSNAME=regular13
EXPNAME=L50I02
F90=pgf90
CDIR=`pwd`
cd ..
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
cp $COMDIR/netlib.f .
cp $COMDIR/common_mtx.f90 .
cp $COMDIR/common_letkf.f90 .
cp $L96DIR/model/lorenz96$ORO.f90 .
cp $L96DIR/obs/h_ope.f90 .
cp $CDIR/letkf.f90 .
$F90 -o letkf SFMT.f90 common.f90 netlib.f common_mtx.f90 common_letkf.f90 lorenz96$ORO.f90 h_ope.f90 letkf.f90
rm *.mod
rm *.o
ln -s $OUTDIR/$OBSNAME/obs.dat .
ln -s $OUTDIR/spinup/init*.dat .
time ./letkf
mkdir -p $OUTDIR/$OBSNAME/$EXPNAME
mv guesmean.dat $OUTDIR/$OBSNAME/$EXPNAME
mv analmean.dat $OUTDIR/$OBSNAME/$EXPNAME
mv gues.dat $OUTDIR/$OBSNAME/$EXPNAME
mv anal.dat $OUTDIR/$OBSNAME/$EXPNAME
cp $CDIR/*.ctl $OUTDIR/$OBSNAME/$EXPNAME

echo "NORMAL END"
