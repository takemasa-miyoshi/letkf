#!/bin/sh
set -e
CDIR=`pwd`
cd ../..
ROMS=`pwd`
DATA1=letkf002/anal/mean
#DATA1=noobs/anal/mean
DATA2=nature
TMPDIR=$CDIR/tmp
#
# Compile
#
rm -rf $TMPDIR
mkdir -p $TMPDIR
cd $TMPDIR
F90=pgf90
PGM=ave_diff
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
ln -s $ROMS/DATA/$DATA1 input1
ln -s $ROMS/DATA/$DATA2 input2
OUTPUT=`echo "$DATA1-$DATA2" | sed -e 's/\//_/g'`
cp input1/0001.nc ave.nc
./$PGM
mv ave.nc $CDIR/ave/$OUTPUT.nc

echo "NORMAL END"
