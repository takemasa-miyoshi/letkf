#!/bin/sh
#=======================================================================
# letkf_cycle.sh
#   To run the SPEEDY-LETKF cycle in parallel computing environment
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
NODE=4
MEMBER=20
OBS=reg3
EXP=M20L500I040
### directory settings
CDIR=`pwd`
cd ../..
SPEEDY=`pwd`
OUTPUT=$SPEEDY/DATA/$OBS/$EXP # data directory
OBSDIR=$SPEEDY/DATA/$OBS/obs  # obs data directory
TMPDIR=$SPEEDY/DATA/tmp/letkf # work directory
LETKF=letkf020.m01
### initial date setting
IYYYY=1982
IMM=01
IDD=01
IHH=00
### final date setting
EYYYY=1982
EMM=05
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
mkdir -p $TMPDIR/ensfcst
cd $TMPDIR/ensfcst
cp $CDIR/ensfcst.sh .
#
# Build SPEEDY model
#
echo '>>>'
echo '>>> BEGIN BUILDING SPEEDY MODEL'
echo '>>>'
cp $SPEEDY/model/source/makefile .
cp $SPEEDY/model/source/*.h .
cp $SPEEDY/model/source/*.f .
cp $SPEEDY/model/source/*.s .

mv par_horres_t30.h atparam.h
mv par_verres.h atparam1.h

cp $SPEEDY/model/ver32.input/cls_*.h .

cp $SPEEDY/model/update/*.h .
cp $SPEEDY/model/update/*.f .
cp $SPEEDY/model/update/makefile .

make imp.exe
echo '>>>'
echo '>>> END BUILDING SPEEDY MODEL'
echo '>>>'
#
# Cycle run ### MAIN LOOP ###
#
while test $IYYYY$IMM$IDD$IHH -le $EYYYY$EMM$EDD$EHH
do
echo '>>>'
echo ">>> BEGIN COMPUTATION OF $IYYYY/$IMM/$IDD/$IHH"
echo '>>>'
TY=`timeinc6hr $IYYYY $IMM $IDD $IHH | cut -c1-4`
TM=`timeinc6hr $IYYYY $IMM $IDD $IHH | cut -c5-6`
TD=`timeinc6hr $IYYYY $IMM $IDD $IHH | cut -c7-8`
TH=`timeinc6hr $IYYYY $IMM $IDD $IHH | cut -c9-10`
#
# LETKF
#
echo " >>"
echo " >> LETKF DATA ASSIMILATION"
echo " >>"
rm -rf $TMPDIR/letkf
mkdir -p $TMPDIR/letkf
cd $TMPDIR/letkf
ln -s $SPEEDY/letkf/$LETKF $LETKF
### inputs
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
ln -s $OUTPUT/gues/$MEM/$IYYYY$IMM$IDD$IHH.grd gs01$MEM.grd
MEM=`expr $MEM + 1`
done
if test -f $OUTPUT/infl_mul/$IYYYY$IMM$IDD$IHH.grd
then
ln -s $OUTPUT/infl_mul/$IYYYY$IMM$IDD$IHH.grd infl_mul.grd
fi
ln -s $OBSDIR/$IYYYY$IMM$IDD$IHH.dat obs01.dat
ln -s $SPEEDY/common/orography_t30.dat fort.21
### mpiexec
mpiexec -n $NODE ./$LETKF < /dev/null
tail -n 17 NOUT-000
### outputs
mv NOUT-000 $OUTPUT/log/$IYYYY$IMM$IDD$IHH.log
if test -f infl_mul.grd
then
cp infl_mul.grd $OUTPUT/infl_mul/$TY$TM$TD$TH.grd
fi
mv gues_me.grd $OUTPUT/gues/mean/$IYYYY$IMM$IDD$IHH.grd
mv gues_sp.grd $OUTPUT/gues/sprd/$IYYYY$IMM$IDD$IHH.grd
mv anal_me.grd $OUTPUT/anal/mean/$IYYYY$IMM$IDD$IHH.grd
mv anal_sp.grd $OUTPUT/anal/sprd/$IYYYY$IMM$IDD$IHH.grd
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
mv anal$MEM.grd $OUTPUT/anal/$MEM/$IYYYY$IMM$IDD$IHH.grd
MEM=`expr $MEM + 1`
done
#
# ensemble forecast
#
echo " >>"
echo " >> ENSEMBLE PREDICTION"
echo " >>"
cd $TMPDIR/ensfcst
MEM=1
while test $MEM -le $MEMBER
do

N=1
while test $N -le $NODE
do
if test $N -lt 10
then
N=0$N
fi

if test $MEM -le $MEMBER
then
if test $MEM -lt 100
then
MEM=0$MEM
fi
if test $MEM -lt 10
then
MEM=0$MEM
fi
echo "MEMBER $MEM in NODE $N"
sh ensfcst.sh $SPEEDY $OUTPUT $IYYYY$IMM$IDD$IHH $TY$TM$TD$TH $MEM $N &
fi
MEM=`expr $MEM + 1`
N=`expr  $N + 1`
done
### wait for the end of parallel processing
time wait
done
#
# Date change ### MAIN LOOP END ###
#
IYYYY=$TY
IMM=$TM
IDD=$TD
IHH=$TH
done

echo "NORMAL END"

