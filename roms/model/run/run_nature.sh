#!/bin/sh
set -e
export OMP_NUM_THREADS=4
CDIR=`pwd`
cd ../..
ROMSDIR=`pwd`
SAVEDIR=$ROMSDIR/DATA/nature
STARTIME=0001
ENDTIME=0720
#INITRST=$ROMSDIR/model/bc/ee6_ini_ee10_rst.0011.nc
INITRST=$ROMSDIR/DATA/spinup/rst04.1.nc
#INITRST=$ROMSDIR/DATA/nature/$STARTIME.nc
WKDIR=$ROMSDIR/tmp/run_nature
rm -rf $WKDIR
mkdir -p $WKDIR
cd $WKDIR
mkdir -p $SAVEDIR
mkdir -p $SAVEDIR/his
### copy
cp $ROMSDIR/model/src/roms .
cp $ROMSDIR/model/bc/*.nc .
cp $ROMSDIR/ncio/chtimestep .
cp $ROMSDIR/model/updates/roms.in.6hr roms.in
### initial
TIME=$STARTIME
cp $INITRST grd.nc
./chtimestep
mv grd.nc $SAVEDIR/$TIME.nc
### cycle run
while test $TIME -lt $ENDTIME
do
echo ">> RUNNING TIME $TIME"
ln -s $SAVEDIR/$TIME.nc input.nc
time ./roms > roms.log
rm input.nc
mv output.0001.nc grd.nc
./chtimestep
TIME=`expr $TIME + 1`
if test $TIME -lt 1000
then
TIME=0$TIME
fi
if test $TIME -lt 100
then
TIME=0$TIME
fi
if test $TIME -lt 10
then
TIME=0$TIME
fi
mv grd.nc $SAVEDIR/$TIME.nc
mv hisout.0000.nc $SAVEDIR/his/$TIME.nc
done

echo "NORMAL END"
