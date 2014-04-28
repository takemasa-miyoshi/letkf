#!/bin/sh
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=16:ompthreads=16
#PBS -o run_nature.log
#PBS -e run_nature.logerr
#PBS -q longp
set -ex
export OMP_NUM_THREADS=16
CDIR=/home/kayoide/enkf/roms/model/run
cd $CDIR
rm -f run_nature.log
rm -f run_nature.logerr
cd ../..
ROMSDIR=`pwd`
SAVEDIR=$ROMSDIR/DATA/nature
IY=2004
IM=01
ID=01
IH=00
EY=2004
EM=03
ED=01
EH=00
source $ROMSDIR/../common/timeinc.sh
#INITRST=$ROMSDIR/model/bc/ee6_ini_ee10_rst.0011.nc
INITRST=$ROMSDIR/DATA/spinup/rst04.1.nc
#INITRST=$ROMSDIR/DATA/nature/$IY$IM$ID$IH.nc
WKDIR=$ROMSDIR/DATA/wkdir/nature
rm -rf $WKDIR
mkdir -p $WKDIR
cd $WKDIR
mkdir -p $SAVEDIR
### copy
cp $ROMSDIR/model/src/roms .
cp $ROMSDIR/model/bc/*.nc .
cp $ROMSDIR/ncio/chtimestep .
cp $ROMSDIR/model/updates/roms.in.6hr roms.in
### initial
cp $INITRST grd.nc
./chtimestep
mv grd.nc $SAVEDIR/$IY$IM$ID${IH}_rst.nc
### cycle run
while test $IY$IM$ID$IH -lt $EY$EM$ED$EH
do
echo ">> RUNNING TIME $IY/$IM/$ID/$IH"
ln -s $SAVEDIR/$IY$IM$ID${IH}_rst.nc input.nc
time ./roms > roms.log
rm input.nc
mv output.0001.nc grd.nc
./chtimestep
### 6hr increment
TY=`timeinc6hr $IY $IM $ID $IH | cut -c1-4`
TM=`timeinc6hr $IY $IM $ID $IH | cut -c5-6`
TD=`timeinc6hr $IY $IM $ID $IH | cut -c7-8`
TH=`timeinc6hr $IY $IM $ID $IH | cut -c9-10`
IY=$TY
IM=$TM
ID=$TD
IH=$TH
### save
mv grd.nc $SAVEDIR/$IY$IM$ID${IH}_rst.nc
mv hisout.0000.nc $SAVEDIR/$IY$IM$ID${IH}_his.nc
done

echo "NORMAL END"
