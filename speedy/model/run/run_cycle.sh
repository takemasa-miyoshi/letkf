#!/bin/sh
#=======================================================================
# run_cycle
#   To run the SPEEDY model for forecast-forecast cycle, useful to
#   generate a nature run.
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
cd ../..
SPEEDY=`pwd`
OUTPUT=$SPEEDY/DATA/nature # directory for output data
TMPDIR=$SPEEDY/model/tmp   # work directory
### initial date setting
IYYYY=1982
IMM=01
IDD=01
IHH=00
### final date setting
EYYYY=1982
EMM=03
EDD=01
EHH=00
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
source $SPEEDY/../common/timeinc.sh
#
# Work directory
#
rm -rf $TMPDIR
mkdir -p $TMPDIR
cd $TMPDIR
#
# Build SPEEDY model
#
echo '>>>BEGIN BUILDING SPEEDY MODEL'
cp $SPEEDY/model/source/makefile .
cp $SPEEDY/model/source/*.h .
cp $SPEEDY/model/source/*.f .
cp $SPEEDY/model/source/*.s .

mv par_horres_t30.h atparam.h
mv par_verres.h atparam1.h

cp $SPEEDY/model/ver32.input/cls_*.h .
cp $SPEEDY/model/ver32.input/inpfiles.s .

cp $SPEEDY/model/update/*.h .
cp $SPEEDY/model/update/*.f .
cp $SPEEDY/model/update/makefile .

make imp.exe

sh inpfiles.s t30

echo '>>>END BUILDING SPEEDY MODEL'
#
# Cycle run ### MAIN LOOP ###
#
while test $IYYYY$IMM$IDD$IHH -le $EYYYY$EMM$EDD$EHH
do
### run
FORT2=1111
ln -fs $OUTPUT/$IYYYY$IMM$IDD$IHH.grd fort.90
echo ">>>BEGIN COMPUTATION OF $IYYYY/$IMM/$IDD/$IHH"
echo $FORT2 > fort.2
echo $IYYYY >> fort.2
echo $IMM >> fort.2
echo $IDD >> fort.2
echo $IHH >> fort.2
time ./imp.exe > out.lis
### date change
TY=`timeinc6hr $IYYYY $IMM $IDD $IHH | cut -c1-4`
TM=`timeinc6hr $IYYYY $IMM $IDD $IHH | cut -c5-6`
TD=`timeinc6hr $IYYYY $IMM $IDD $IHH | cut -c7-8`
TH=`timeinc6hr $IYYYY $IMM $IDD $IHH | cut -c9-10`
IYYYY=$TY
IMM=$TM
IDD=$TD
IHH=$TH
### output
mv $IYYYY$IMM$IDD${IHH}*.grd $OUTPUT
done

echo "NORMAL END"

