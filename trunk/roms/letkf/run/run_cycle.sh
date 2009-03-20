#!/bin/sh
#=======================================================================
# letkf_cycle.sh
#   To run the ROMS-LETKF cycle in parallel computing environment
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
NODE=4
MEMBER=20
### directory settings
CDIR=`pwd`
cd ../..
ROMS=`pwd`
OUTPUT=$ROMS/DATA/letkf001  # data directory
OBSDIR=$ROMS/DATA/obs   # obs data directory
TMPDIR=$ROMS/letkf/tmp  # work directory
LETKF=letkf020.m01
### initial date setting
IT=0001
### final date setting
ET=0120
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
while test $IT -le $ET
do
echo '>>>'
echo ">>> BEGIN COMPUTATION OF $IT"
echo '>>>'
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
ln -s $OUTPUT/gues/$MEM/$IT.nc gs01$MEM.nc
cp $OUTPUT/gues/$MEM/$IT.nc anal$MEM.nc
MEM=`expr $MEM + 1`
done
ln -s $OBSDIR/$IT.dat obs01.dat
ln -s $ROMS/model/bc/ee6_grd.nc grd.nc
cp anal001.nc gues_me.nc
cp anal001.nc gues_sp.nc
cp anal001.nc anal_me.nc
cp anal001.nc anal_sp.nc
### mpiexec
mpiexec -n $NODE ./$LETKF < /dev/null
tail -n 17 NOUT-000
cp NOUT-000 $OUTPUT/log/$IT.log
### outputs
mv gues_me.nc $OUTPUT/gues/mean/$IT.nc
mv gues_sp.nc $OUTPUT/gues/sprd/$IT.nc
mv anal_me.nc $OUTPUT/anal/mean/$IT.nc
mv anal_sp.nc $OUTPUT/anal/sprd/$IT.nc
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
mv anal$MEM.nc $OUTPUT/anal/$MEM/$IT.nc
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
sh ensfcst.sh $ROMS $OUTPUT $IT $MEM $N &
fi
MEM=`expr $MEM + 1`
N=`expr  $N + 1`
done
### wait for the end of parallel processing
time wait
done
#
# Time change ### MAIN LOOP END ###
#
IT=`expr $IT + 1`
if test $IT -lt 1000
then
IT=0$IT
fi
if test $IT -lt 100
then
IT=0$IT
fi
if test $IT -lt 10
then
IT=0$IT
fi
done

echo "NORMAL END"

