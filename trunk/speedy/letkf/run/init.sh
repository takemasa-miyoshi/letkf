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
cd ../..
SPEEDY=`pwd`
NATURE=$SPEEDY/DATA/nature # nature run
OUTPUT=$SPEEDY/DATA/test   # directory for new experiment
### initial date setting
IYYYY=1982
IMM=01
IDD=01
IHH=00
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
source $SPEEDY/../common/timeinc.sh
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
mkdir -p $OUTPUT/anal_f/$MEM
mkdir -p $OUTPUT/gues/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/anal/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/anal_f/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/gues/$MEM
MEM=`expr $MEM + 1`
done

for MEM in mean sprd
do
mkdir -p $OUTPUT/anal/$MEM
mkdir -p $OUTPUT/anal_f/$MEM
mkdir -p $OUTPUT/gues/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/anal/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/anal_f/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/gues/$MEM
done
### copy initial conditions
TY=1982
TM=02
TD=01
TH=00
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
cp $NATURE/$TY$TM$TD$TH.grd $OUTPUT/gues/$MEM/$IYYYY$IMM$IDD$IHH.grd
UY=`timeinc6hr $TY $TM $TD $TH | cut -c1-4`
UM=`timeinc6hr $TY $TM $TD $TH | cut -c5-6`
UD=`timeinc6hr $TY $TM $TD $TH | cut -c7-8`
UH=`timeinc6hr $TY $TM $TD $TH | cut -c9-10`
TY=`timeinc6hr $UY $UM $UD $UH | cut -c1-4`
TM=`timeinc6hr $UY $UM $UD $UH | cut -c5-6`
TD=`timeinc6hr $UY $UM $UD $UH | cut -c7-8`
TH=`timeinc6hr $UY $UM $UD $UH | cut -c9-10`
MEM=`expr $MEM + 1`
done

echo "NORMAL END"
