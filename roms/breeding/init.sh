#!/bin/sh
#=======================================================================
# init.sh
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
MEMBER=3
### directory settings
CDIR=`pwd`
cd ..
ROMS=`pwd`
SPINUP=$ROMS/DATA/spinup
NATURE=$ROMS/DATA/nature
OUTPUT=$ROMS/DATA/breeding   # directory for new experiment
### initial time setting
IT=0001
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
cd $CDIR
### clean
rm -rf $OUTPUT
### mkdir
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
mkdir -p $OUTPUT/fcst/$MEM
mkdir -p $OUTPUT/bv/$MEM
MEM=`expr $MEM + 1`
done
### copy initial conditions
I=1
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
cp $SPINUP/rst03.$I.nc $OUTPUT/fcst/$MEM/$IT.nc
ln -s $NATURE/$IT.nc in.nc
ln -s $OUTPUT/fcst/$MEM/$IT.nc out.nc
$ROMS/ncio/choceantime
rm in.nc
rm out.nc
I=`expr $I + 1`
MEM=`expr $MEM + 1`
done

echo "NORMAL END"
