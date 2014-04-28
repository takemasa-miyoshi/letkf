#!/bin/sh
set -e
CDIR=`pwd`
cd ..
ROMS=`pwd`
DATA=$ROMS/DATA/spinup
OUTPUT=spup_uv.dat
START=0
END=4
#
# Compile
#
cd $CDIR
F90=pgf90
PGM=diag_spup_uv
LIB_NETCDF="-L/usr/local/netcdf/lib -lnetcdf"
INC_NETCDF="-I/usr/local/netcdf/include"
ln -fs ../../common/SFMT.f90 .
ln -fs ../../common/common.f90 .
ln -fs ../common/common_roms.f90 .
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
rm -f $OUTPUT
IY=$START
while test $IY -le $END
do
echo ">> PROCESSING YEAR 0$IY"
for DATE in 0000 0030 0060 0090 0120 0150 0180 0210 0240 0270 0300 0330
do
ln -s $DATA/hist0$IY.$DATE.nc input.nc
./$PGM >> $OUTPUT
rm input.nc
done
IY=`expr $IY + 1`
done
rm $PGM

echo "NORMAL END"
