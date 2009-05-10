#!/bin/sh
set -e
CDIR=`pwd`
cd ../..
ROMS=`pwd`
#DATA="nature letkf001/anal/sprd letkf001/gues/sprd bv_ubar_5E-3/bv/001"
DATA="letkf002/anal/sprd"
TMPDIR=$CDIR/tmp
#
# Compile
#
rm -rf $TMPDIR
mkdir -p $TMPDIR
cd $TMPDIR
F90=pgf90
PGM=ave
LIB_NETCDF="-L/usr/local/netcdf/lib -lnetcdf"
INC_NETCDF="-I/usr/local/netcdf/include"
ln -fs $ROMS/../common/SFMT.f90 .
ln -fs $ROMS/../common/common.f90 .
ln -fs $ROMS/common/common_roms.f90 .
ln -fs $CDIR/$PGM.f90 .
F90OPT='-byteswapio -tp p7-64 -fast -O3'
$F90 $F90OPT -c SFMT.f90
$F90 $F90OPT -c common.f90
$F90 $F90OPT $INC_NETCDF -c common_roms.f90
$F90 $F90OPT $INC_NETCDF -c $PGM.f90
$F90 -o $PGM *.o $LIB_NETCDF
rm SFMT.f90
rm common.f90
rm common_roms.f90
rm *.mod
rm *.o
#
# main loop
#
for DATADIR in $DATA
do
ln -s $ROMS/DATA/$DATADIR input
OUTPUT=`echo $DATADIR | sed -e 's/\//_/g'`
cp input/0001.nc ave.nc
./$PGM
rm input
mv ave.nc $CDIR/ave/$OUTPUT.nc
done

echo "NORMAL END"
