#!/bin/sh
#=======================================================================
# init.sh
#   This script prepares for new LETKF cycle-run experiment
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
MEMBER=20
### directory settings
CDIR=`pwd`
cd ../..
ROMS=`pwd`
OBS=5x20_new
EXP=LETKF20H50V050
NATURE=$ROMS/DATA/nature # nature run
SPINUP=$ROMS/DATA/spinup # spin-up run
OUTPUT=$ROMS/DATA/$OBS/$EXP   # directory for new experiment
### initial time setting
IT=2004010100
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
cd $CDIR
### clean
rm -rf $OUTPUT
### mkdir
mkdir -p $OUTPUT/log
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
mkdir -p $OUTPUT/anal/$MEM
mkdir -p $OUTPUT/gues/$MEM
MEM=`expr $MEM + 1`
done

for MEM in mean sprd
do
mkdir -p $OUTPUT/anal/$MEM
mkdir -p $OUTPUT/gues/$MEM
done

mkdir -p $OUTPUT/infl_mul
### copy initial conditions
Y=0
MEM=0
while test $Y -le 3
do
I=1
if test $Y -eq 0
then
I=4
fi
while test $I -le 6
do
MEM=`expr $MEM + 1`
if test $MEM -gt $MEMBER
then
break
fi
if test $MEM -lt 100
then
MEM=0$MEM
fi
if test $MEM -lt 10
then
MEM=0$MEM
fi
echo "copying.. rst0$Y.$I.nc --> member $MEM"
cp $SPINUP/rst0$Y.$I.nc $OUTPUT/gues/$MEM/${IT}_rst.nc
ln -s $NATURE/${IT}_rst.nc in.nc
ln -s $OUTPUT/gues/$MEM/${IT}_rst.nc out.nc
$ROMS/ncio/choceantime
rm in.nc
rm out.nc
I=`expr $I + 1`
done
if test $MEM -gt $MEMBER
then
break
fi
Y=`expr $Y + 1`
done

cp $NATURE/../infl_mul.nc $OUTPUT/infl_mul/${IT}_rst.nc

echo "NORMAL END"
