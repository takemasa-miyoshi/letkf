#!/bin/sh
set -e
CDIR=`pwd`
cd ..
ROMS=`pwd`
export OMP_NUM_THREADS=4
INC=1
#1:6hr, 4:1d, 20:5d
NATURE=$ROMS/DATA/nature
SAVEDIR=$ROMS/DATA/breeding
#SAVEDIR=$ROMS/DATA/bv_ubar_1E-2
WKDIR=$CDIR/tmp
MEMBER=1
START=0001
END=0240
#END=0720
#
# Compile
#
cd $CDIR
F90=pgf90
PGM=breeding
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
# main process
#
rm -rf $WKDIR
mkdir -p $WKDIR
cd $WKDIR
cp $ROMS/model/src/roms .
cp $ROMS/model/bc/*.nc .
cp $ROMS/ncio/chtimestep .
if test $INC -eq 1
then
cp $ROMS/model/updates/roms.in.6hr roms.in
elif test $INC -eq 4
then
cp $ROMS/model/updates/roms.in.1dy roms.in
elif test $INC -eq 20
then
cp $ROMS/model/updates/roms.in.5dy roms.in
else
echo "ERROR in INC=$INC"
exit
fi
mv $CDIR/$PGM .
#
# main loop
#
TIME=$START
while test $TIME -le $END
do
echo ">> RUNNING TIME $TIME"
TIME2=`expr $TIME + $INC`
if test $TIME2 -lt 1000
then
TIME2=0$TIME2
fi
if test $TIME2 -lt 100
then
TIME2=0$TIME2
fi
if test $TIME2 -lt 10
then
TIME2=0$TIME2
fi
MEM=1
while test $MEM -le $MEMBER
do
if test $MEM -lt 100
then
MEM=0$MEM
fi
if test $MEM -lt 10
then
MEM=0$MEM
fi
ln -s $NATURE/$TIME.nc nature.nc       #(in)    nature run
cp $SAVEDIR/fcst/$MEM/$TIME.nc ptb.nc  #(inout) perturbed fcst -> ptb ic
cp ptb.nc bv.nc                        #(out)   bv (before normalization)
./$PGM >> $SAVEDIR/rescale$MEM.dat
rm nature.nc
mv bv.nc $SAVEDIR/bv/$MEM/$TIME.nc
mv ptb.nc grd.nc
./chtimestep
mv grd.nc input.nc
time ./roms > roms.log
rm input.nc
mv output.0001.nc $SAVEDIR/fcst/$MEM/$TIME2.nc
MEM=`expr $MEM + 1`
done
TIME=$TIME2
done

echo "NORMAL END"
