#!/bin/sh
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=16:ompthreads=16
#PBS -o run_spinup.log
#PBS -e run_spinup.logerr
#PBS -q longp
set -ex
export OMP_NUM_THREADS=16
CDIR=/home/kayoide/enkf/roms/model/run
cd $CDIR
rm -f run_spinup.log
rm -f run_spinup.logerr
cd ../..
ROMSDIR=`pwd`
SAVEDIR=$ROMSDIR/DATA/spinup
STARTYEAR=00
ENDYEAR=05
INITRST=$ROMSDIR/model/bc/ee6_ini_ee10_rst.0011.nc
#WKDIR=/workp/bco/kayoide/tmp/run_spinup
#WKDIR=/lscratch/kayoide/run_spinup
WKDIR=$ROMSDIR/DATA/wkdir/spinup
rm -rf $WKDIR
mkdir -p $WKDIR
cd $WKDIR
mkdir -p $SAVEDIR
### copy
cp $ROMSDIR/model/src/roms .
cp $ROMSDIR/model/bc/*.nc .
cp $ROMSDIR/ncio/chtimestep .
cp $ROMSDIR/model/updates/roms.in.1yr roms.in
### initial
YEAR=$STARTYEAR
cp $INITRST grd.nc
./chtimestep
mv grd.nc $SAVEDIR/rst$YEAR.1.nc
### cycle run
while test $YEAR -lt $ENDYEAR
do
echo ">> START RUNNING YEAR $YEAR"
ln -s $SAVEDIR/rst$YEAR.1.nc input.nc
time ./roms > roms.log
rm input.nc
I=1
while test $I -lt 6
do
mv output.000$I.nc grd.nc
./chtimestep
I=`expr $I + 1`
mv grd.nc $SAVEDIR/rst$YEAR.$I.nc
done
for HISTDAY in 0000 0030 0060 0090 0120 0150 0180 0210 0240 0270 0300 0330
do
mv history.$HISTDAY.nc $SAVEDIR/hist$YEAR.$HISTDAY.nc
done
mv output.0006.nc grd.nc
./chtimestep
YEAR=`expr $YEAR + 1`
if test $YEAR -lt 10
then
YEAR=0$YEAR
fi
mv grd.nc $SAVEDIR/rst$YEAR.1.nc
done

mv roms.log $SAVEDIR
rm -rf $WKDIR
echo "NORMAL END"
