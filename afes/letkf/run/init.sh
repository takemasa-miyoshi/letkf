#!/bin/sh
#=======================================================================
# init.sh
#   This script prepares for new LETKF cycle-run experiment
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
MEMBER=10
OBS=test
EXP=M10
### directory settings
cd ../..
AFES=`pwd`
OUTPUT=/S/data04/G4015/y0266/alera2/$OBS/$EXP  # directory for new experiment
INDIR=/S/data04/G4015/y0260/alera2/sample
### initial date setting
IYYYY=1982
IMM=01
IDD=01
IHH=00
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
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
mkdir -p $OUTPUT/anal/$MEM
mkdir -p $OUTPUT/gues/$MEM
cp $AFES/common/yyyymmddhh.ctl $OUTPUT/anal/$MEM
cp $AFES/common/yyyymmddhh.ctl $OUTPUT/gues/$MEM
MEM=`expr $MEM + 1`
done

for MEM in mean sprd
do
mkdir -p $OUTPUT/anal/$MEM
mkdir -p $OUTPUT/gues/$MEM
cp $AFES/common/yyyymmddhh.ctl $OUTPUT/anal/$MEM
cp $AFES/common/yyyymmddhh.ctl $OUTPUT/gues/$MEM
done
### copy initial conditions
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
cp $INDIR/0$MEM/IC.grd $OUTPUT/gues/$MEM/$IYYYY$IMM$IDD$IHH.grd
MEM=`expr $MEM + 1`
done

echo "NORMAL END"
