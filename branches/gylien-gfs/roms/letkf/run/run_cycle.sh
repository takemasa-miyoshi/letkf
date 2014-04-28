#!/bin/sh
#PBS -l walltime=12:00:00
#PBS -l ncpus=40
#PBS -o run_cycle.log
#PBS -e run_cycle.logerr
#PBS -q longp
#=======================================================================
# letkf_cycle.sh
#   To run the ROMS-LETKF cycle in parallel computing environment
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
ANODE=40
FNODE=20
export OMP_NUM_THREADS=2
MEMBER=20
### directory settings
CDIR=/home/kayoide/enkf/roms/letkf/run
cd $CDIR
rm -f run_cycle.log
rm -f run_cycle.logerr
cd ../..
ROMS=`pwd`
source $ROMS/../common/timeinc.sh
OBS=5x20_new
EXP=LETKF20H50V050
if test $MEMBER -lt 100
then
MEM=0$MEMBER
else
MEM=$MEMBER
fi
OUTPUT=$ROMS/DATA/$OBS/$EXP  # data directory
OBSDIR=$ROMS/DATA/$OBS/letkf_obs   # obs data directory
#TMPDIR=/lscratch/kayoide/letkf  # work directory
TMPDIR=$ROMS/DATA/wkdir/$EXP # work directory
LETKF=letkf020.m01
### initial date setting
IY=2004
IM=01
ID=01
IH=00
### final date setting
EY=2004
EM=02
ED=01
EH=00
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
#
# Work directory
#
rm -rf $TMPDIR
mkdir -p $TMPDIR/ensfcst
cd $TMPDIR/ensfcst
cp $CDIR/ensfcst.sh .
#
# Cycle run ### MAIN LOOP ###
#
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do
echo '>>>'
echo ">>> BEGIN COMPUTATION OF $IY/$IM/$ID/$IH"
echo '>>>'
TY=`timeinc6hr $IY $IM $ID $IH | cut -c1-4`
TM=`timeinc6hr $IY $IM $ID $IH | cut -c5-6`
TD=`timeinc6hr $IY $IM $ID $IH | cut -c7-8`
TH=`timeinc6hr $IY $IM $ID $IH | cut -c9-10`
#
# LETKF
#
echo " >>"
echo " >> LETKF DATA ASSIMILATION"
echo " >>"
rm -rf $TMPDIR/letkf
mkdir -p $TMPDIR/letkf
cd $TMPDIR/letkf
ln -s $ROMS/letkf/$LETKF $LETKF
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
ln -s $OUTPUT/gues/$MEM/$IY$IM$ID${IH}_rst.nc gs01$MEM.nc
cp $OUTPUT/gues/$MEM/$IY$IM$ID${IH}_rst.nc anal$MEM.nc
MEM=`expr $MEM + 1`
done
#for OBTYPE in hfradar_uv ship_sst ship_ssh profiles prof_uv
for OBTYPE in hfradar_uv ship_sst profiles
do
cat $OBSDIR/$OBTYPE/$IY$IM$ID${IH}_obs.dat >> obs01.dat
done
ln -s $ROMS/model/bc/ee6_grd.nc grd.nc
cp anal001.nc gues_me.nc
cp anal001.nc gues_sp.nc
cp anal001.nc anal_me.nc
cp anal001.nc anal_sp.nc
if test -f $OUTPUT/infl_mul/$IY$IM$ID${IH}_rst.nc
then
cp $OUTPUT/infl_mul/$IY$IM$ID${IH}_rst.nc infl_mul.nc
fi
### mpiexec
mpirun -np $ANODE ./$LETKF < /dev/null
tail -n 17 NOUT-000
cp NOUT-000 $OUTPUT/log/$IY$IM$ID${IH}.log
### outputs
if test -f infl_mul.nc
then
cp infl_mul.nc $OUTPUT/infl_mul/$TY$TM$TD${TH}_rst.nc
fi
mv gues_me.nc $OUTPUT/gues/mean/$IY$IM$ID${IH}_rst.nc
mv gues_sp.nc $OUTPUT/gues/sprd/$IY$IM$ID${IH}_rst.nc
mv anal_me.nc $OUTPUT/anal/mean/$IY$IM$ID${IH}_rst.nc
mv anal_sp.nc $OUTPUT/anal/sprd/$IY$IM$ID${IH}_rst.nc
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
mv anal$MEM.nc $OUTPUT/anal/$MEM/$IY$IM$ID${IH}_rst.nc
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
while test $N -le $FNODE
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
sh ensfcst.sh $ROMS $OUTPUT $IY$IM$ID$IH $TY$TM$TD$TH $MEM $N &
fi
MEM=`expr $MEM + 1`
N=`expr  $N + 1`
done
### wait for the end of parallel processing
time wait
done
#
# Remove data
#
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
rm $OUTPUT/gues/$MEM/$IY$IM$ID${IH}_rst.nc
rm $OUTPUT/anal/$MEM/$IY$IM$ID${IH}_rst.nc
MEM=`expr $MEM + 1`
done
#
# Time change ### MAIN LOOP END ###
#
IY=$TY
IM=$TM
ID=$TD
IH=$TH
done

echo "NORMAL END"

