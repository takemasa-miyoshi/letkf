#!/bin/sh
#=======================================================================
# run_first
#   To run the SPEEDY model for the first time. That is, you do not have
#   a gridded initial condition, so that the model starts from the
#   atmosphere at rest. The atmosphere at rest means zero winds
#   everywhere and constant T with vertical profile of the standard
#   atmosphere.
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
cd ../..
SPEEDY=`pwd`
OUTPUT=$SPEEDY/DATA/nature # directory for output data
TMPDIR=$SPEEDY/model/tmp   # work directory
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
mkdir -p $OUTPUT
cp $SPEEDY/common/yyyymmddhh.ctl $OUTPUT
cp $SPEEDY/common/yyyymmddhh_p.ctl $OUTPUT
#
# Work directory
#
rm -rf $TMPDIR
mkdir -p $TMPDIR
cd $TMPDIR
#
# 1-year run
#
### date setting
IYYYY=1981
IMM=01
IDD=01
IHH=00
### build
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
cp $SPEEDY/model/update/cls_instep.h_1yr cls_instep.h # for 1-yr run

make imp.exe

sh inpfiles.s t30

echo '>>>END BUILDING SPEEDY MODEL'
### run
FORT2=0
echo ">>>BEGIN 1-YEAR COMPUTATION OF $IYYYY/$IMM/$IDD/$IHH"
echo $FORT2 > fort.2
echo $IYYYY >> fort.2
echo $IMM >> fort.2
echo $IDD >> fort.2
echo $IHH >> fort.2
time ./imp.exe > out.lis
mv fort.10 restart.dat
### clean up
rm -f *.o
rm -f imp.exe
#
# 6-hr run
#
### date setting
IYYYY=1981
IMM=12
IDD=31
IHH=00
### build
echo '>>>BEGIN BUILDING SPEEDY MODEL'
cp $SPEEDY/model/update/cls_instep.h ./ # for 6-hr run
make imp.exe
### run 1
FORT2=1110
ln -s restart.dat fort.3
echo ">>>BEGIN COMPUTATION OF $IYYYY/$IMM/$IDD/$IHH"
echo $FORT2 > fort.2
echo $IYYYY >> fort.2
echo $IMM >> fort.2
echo $IDD >> fort.2
echo $IHH >> fort.2
time ./imp.exe > out.lis
IHH=06
FORT2=1111
mv $IYYYY$IMM$IDD$IHH.grd $OUTPUT
mv $IYYYY$IMM$IDD${IHH}_p.grd $OUTPUT
### run 2
ln -fs $OUTPUT/$IYYYY$IMM$IDD$IHH.grd fort.90
echo ">>>BEGIN COMPUTATION OF $IYYYY/$IMM/$IDD/$IHH"
echo $FORT2 > fort.2
echo $IYYYY >> fort.2
echo $IMM >> fort.2
echo $IDD >> fort.2
echo $IHH >> fort.2
time ./imp.exe > out.lis
IHH=12
mv $IYYYY$IMM$IDD$IHH.grd $OUTPUT
mv $IYYYY$IMM$IDD${IHH}_p.grd $OUTPUT
### run 3
ln -fs $OUTPUT/$IYYYY$IMM$IDD$IHH.grd fort.90
echo ">>>BEGIN COMPUTATION OF $IYYYY/$IMM/$IDD/$IHH"
echo $FORT2 > fort.2
echo $IYYYY >> fort.2
echo $IMM >> fort.2
echo $IDD >> fort.2
echo $IHH >> fort.2
time ./imp.exe > out.lis
IHH=18
mv $IYYYY$IMM$IDD$IHH.grd $OUTPUT
mv $IYYYY$IMM$IDD${IHH}_p.grd $OUTPUT
### run 4
ln -fs $OUTPUT/$IYYYY$IMM$IDD$IHH.grd fort.90
echo ">>>BEGIN COMPUTATION OF $IYYYY/$IMM/$IDD/$IHH"
echo $FORT2 > fort.2
echo $IYYYY >> fort.2
echo $IMM >> fort.2
echo $IDD >> fort.2
echo $IHH >> fort.2
time ./imp.exe > out.lis
IYYYY=1982
IMM=01
IDD=01
IHH=00
mv $IYYYY$IMM$IDD$IHH.grd $OUTPUT
mv $IYYYY$IMM$IDD${IHH}_p.grd $OUTPUT

echo "NORMAL END"

