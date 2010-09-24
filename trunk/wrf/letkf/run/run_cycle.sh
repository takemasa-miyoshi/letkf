#!/bin/sh
#=======================================================================
# letkf_cycle.sh
#   To run the WRF-LETKF cycle in parallel computing environment
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
MEMBER=20
RSH=rsh
RCP="rcp -q"
### parallel computation setting
NODES="node2 node3 node4 node5 node6 node7 node8" # names of each node
WRF_MEMBER_PER_NODE=3 # members per node for WRF
WRF_PROC=2            # processes per member for WRF
NCPU_PER_NODE=6       # number of CPUs per node for LETKF
### LETKF setting
WINDOW=360            # assimilation window length (min)
GUESFT=540            # first guess forecast length (min)
OBS="prepbufr/ucar"   # "prepbufr/ucar" or "none"
EXP=EXP001            # name of experiment
### initial date setting
IY=2008
IM=09
ID=03
IH=12
IMN=00
### final date setting
EY=2008
EM=09
ED=03
EH=12
EMN=00
### adaptive inflation
ADAPTINFL=1 #1:ON, 0:OFF
### restart setting
ITER=1 #first iteration number, usually 1, but > 1 for restart
### job setting
ENS=1
ANL=1
SUCC_INFL_MUL=
### directory settings
CDIR=`pwd`
cd ../..
WRF=`pwd`
OBSDIR=$WRF/DATA/obs/$OBS          # obs data directory
FNLDIR=$WRF/DATA/FNL               # NCEP/NCAR data directory
TMPDIR=/local/data/$USER/work/$EXP # work directory
LETKF=$WRF/letkf/letkf020.m01      # LETKF module
INITDAT=$WRF/DATA/init_spinup      # initial conditions
NCIO=$WRF/ncio
MPIBIN=/usr/local/mpich2/bin
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
echo '>>>'
echo ">>> STARTING WRF-LETKF CYCLE FROM $IY/$IM/$ID/$IH TO $EY/$EM/$ED/$EH"
echo '>>>'
NOOBS=0
if test $OBS = "none"
then
ADAPTINFL=0
NOOBS=1
fi
MEMBERP=`expr $MEMBER + 1`
NODE0=`echo $NODES | awk '{print $1}'`
NUM_NODE=0
for NODE in $NODES
do
NUM_NODE=`expr $NUM_NODE + 1`
done
LETKF_PROC=`expr $NUM_NODE \* $NCPU_PER_NODE`
cd $CDIR
source util.sh
#
# Work directory
#
if test $ITER -eq 1
then
echo "Setting up work directories.."
echo " >> removing old work directories.."
rm -rf $TMPDIR
mkdir -p $TMPDIR
for NODE in $NODES
do
$RSH $NODE rm -rf $TMPDIR &
done
time wait
for NODE in $NODES
do
$RSH $NODE mkdir -p $TMPDIR/gues
$RSH $NODE mkdir -p $TMPDIR/anal
$RSH $NODE mkdir -p $TMPDIR/model
$RSH $NODE cp -r $NCIO $TMPDIR
$RSH $NODE cp $CDIR/ensfcst.sh $TMPDIR/model
$RSH $NODE cp $WRF/letkf/run/util.sh $TMPDIR/model
$RSH $NODE cp -RL $WRF/model/WPS     $TMPDIR/model &
$RSH $NODE cp -RL $WRF/model/WRF     $TMPDIR/model &
$RSH $NODE cp -RL $WRF/model/ARWpost $TMPDIR/model &
$RSH $NODE cp -RL $WRF/model/GEOG    $TMPDIR/model &
done
echo " >> copying WRF files.."
time wait
fi
# 
# Cycle run ### MAIN LOOP ###
#
ITER=1
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do

echo '>>>'
echo ">>> BEGIN COMPUTATION OF $IY/$IM/$ID/$IH"
echo '>>>'

for NODE in $NODES
do
$RCP $CDIR/ensfcst.sh $NODE:$TMPDIR/model
$RCP $CDIR/util.sh    $NODE:$TMPDIR/model
done

cd $TMPDIR
date_edit $IY $IM $ID $IH 00 $GUESFT > fcst_time.txt
read FY FM FD FH FMN < fcst_time.txt
rm -f fcst_time.txt

date_edit $IY $IM $ID $IH 00 $WINDOW > anal_time.txt
read AY AM AD AH AMN < anal_time.txt
rm -f anal_time.txt
#
# ensemble forecast
#
if test $ENS -eq 1
then
echo " >>"
echo " >> ENSEMBLE FORECASTS"
echo " >>"
### mpdboot
for NODE in $NODES
do
$RSH $NODE $MPIBIN/mpdboot -n 1
done
###
ENDMEMBER=0
M=1
while test $ENDMEMBER -eq 0
do
ILOOP=1
while test $ILOOP -le $WRF_MEMBER_PER_NODE -a $ENDMEMBER -eq 0
do
for NODE in $NODES
do
if test $M -lt 10
then
MEM=00$M
else
MEM=0$M
fi
echo "  > FORECAST MEMBER $MEM in NODE $NODE"
if test $ITER -eq 1
then
FLAG=2 # use forecast data as initial file
ANALDAT=$INITDAT/$MEM/wrfout_d01_${IY}-${IM}-${ID}_${IH}:${IMN}:00
else
FLAG=3 #3: based on NCEP Reanalysis, 4: based on previous forecast
ANALDAT=$TMPDIR/anal/${IY}${IM}${ID}${IH}${IMN}/anal${MEM}.grd
fi
if test $M -eq $MEMBERP
then
ANALDAT=$TMPDIR/anal/${IY}${IM}${ID}${IH}${IMN}/anal_me.grd
ENDMEMBER=1
fi
$RSH $NODE sh $TMPDIR/model/ensfcst.sh $TMPDIR $IY$IM$ID$IH $FY$FM$FD$FH $MEM $WRF_PROC $FNLDIR $ANALDAT $FLAG &
M=`expr $M + 1`
if test $ENDMEMBER -eq 1
then
break
fi
done
ILOOP=`expr $ILOOP + 1`
done
### wait for parallel processing
time wait
done
### mpdallexit
for NODE in $NODES
do
$RSH $NODE $MPIBIN/mpdallexit
done
fi
#
# LETKF
#
if test $ANL -eq 1
then
echo " >>"
echo " >> LETKF DATA ASSIMILATION"
echo " >>"
# set up wkdir
for NODE in $NODES
do
$RSH $NODE rm -rf $TMPDIR/letkf
$RSH $NODE mkdir -p $TMPDIR/letkf
$RSH $NODE cp $LETKF $TMPDIR/letkf
if test $NOOBS -eq 1
then
$RSH $NODE touch $TMPDIR/letkf/obs01.dat
$RSH $NODE touch $TMPDIR/letkf/obs02.dat
$RSH $NODE touch $TMPDIR/letkf/obs03.dat
$RSH $NODE touch $TMPDIR/letkf/obs04.dat
$RSH $NODE touch $TMPDIR/letkf/obs05.dat
$RSH $NODE touch $TMPDIR/letkf/obs06.dat
$RSH $NODE touch $TMPDIR/letkf/obs07.dat
else
$RCP $OBSDIR/obs$AY$AM$AD$AH/t-3.dat $NODE:$TMPDIR/letkf/obs01.dat
$RCP $OBSDIR/obs$AY$AM$AD$AH/t-2.dat $NODE:$TMPDIR/letkf/obs02.dat
$RCP $OBSDIR/obs$AY$AM$AD$AH/t-1.dat $NODE:$TMPDIR/letkf/obs03.dat
$RCP $OBSDIR/obs$AY$AM$AD$AH/t.dat   $NODE:$TMPDIR/letkf/obs04.dat
$RCP $OBSDIR/obs$AY$AM$AD$AH/t+1.dat $NODE:$TMPDIR/letkf/obs05.dat
$RCP $OBSDIR/obs$AY$AM$AD$AH/t+2.dat $NODE:$TMPDIR/letkf/obs06.dat
$RCP $OBSDIR/obs$AY$AM$AD$AH/t+3.dat $NODE:$TMPDIR/letkf/obs07.dat
fi
done
# infl_mul (1 node)
if test $ITER -gt 1 -a $ADAPTINFL -eq 1
then
echo "  > Copying infl_mul.grd"
$RSH $NODE0 cp $TMPDIR/anal/$IY$IM$ID$IH$IMN/infl_mul.grd $TMPDIR/letkf
elif test $ITER -eq 1 -a "X$SUCC_INFL_MUL" != "X"
then
echo "  > Copying $SUCC_INFL_MUL"
$RSH $NODE0 cp $SUCC_INFL_MUL $TMPDIR/letkf
fi
# first guess
ENDMEMBER=0
M=1
while test $ENDMEMBER -eq 0
do
for NODE in $NODES
do
if test $M -eq $MEMBERP
then
NODE1=$NODE
ENDMEMBER=1
break
fi
if test $M -lt 10
then
MEM=00$M
else
MEM=0$M	
fi 
$RSH $NODE ln -s  $TMPDIR/gues/$MEM/$AY$AM$AD$AH$AMN/gs* $TMPDIR/letkf
$RSH $NODE ln -fs  $TMPDIR/gues/$MEM/$AY$AM$AD$AH$AMN/const.grd $TMPDIR/letkf
M=`expr $M + 1`
done
done
# make a machine file for MPICH2
rm -f mpd.hosts
for NODE in $NODES
do
echo $NODE >> mpd.hosts
done
$RCP mpd.hosts $NODE0:$TMPDIR/letkf
M=1
rm -f machinefile
while test $M -le $LETKF_PROC
do
for NODE in $NODES
do
echo $NODE >> machinefile
M=`expr $M + 1`
done
done
$RCP machinefile $NODE0:$TMPDIR/letkf
# make a script for LETKF
cat << EOF > letkf.sh
#!/bin/sh
cd $TMPDIR/letkf
$MPIBIN/mpdboot -n $NUM_NODE -f ./mpd.hosts
$MPIBIN/mpdtrace
$MPIBIN/mpiexec -machinefile ./machinefile -n $LETKF_PROC ./`basename $LETKF` < /dev/null
$MPIBIN/mpdallexit
EOF

echo "  > Starting LETKF"
$RCP letkf.sh $NODE0:$TMPDIR/letkf
$RSH $NODE0 sh $TMPDIR/letkf/letkf.sh
echo "  > END of LETKF"

for NODE in $NODES
do
$RSH $NODE mkdir -p $TMPDIR/anal/$AY$AM$AD$AH$AMN
$RSH $NODE mv $TMPDIR/letkf/anal*.grd $TMPDIR/anal/$AY$AM$AD$AH$AMN
$RSH $NODE mv $TMPDIR/letkf/NOUT-*    $TMPDIR/anal/$AY$AM$AD$AH$AMN
if test $NODE = $NODE0
then
$RSH $NODE mv $TMPDIR/letkf/gues_??.grd  $TMPDIR/anal/$AY$AM$AD$AH$AMN
if test $ADAPTINFL -eq 1
then
$RSH $NODE mv $TMPDIR/letkf/infl_mul.grd $TMPDIR/anal/$AY$AM$AD$AH$AMN
fi
if test $NOOBS -eq 0
then
$RSH $NODE mv $TMPDIR/letkf/obs.dat $TMPDIR/anal/$AY$AM$AD$AH$AMN
fi
$RSH $NODE1 mkdir -p $TMPDIR/anal/$AY$AM$AD$AH$AMN
$RSH $NODE $RCP $TMPDIR/anal/$AY$AM$AD$AH$AMN/anal_me.grd $NODE1:$TMPDIR/anal/$AY$AM$AD$AH$AMN
fi
done

fi
#
# Time change ### MAIN LOOP END ###
#
IY=$AY
IM=$AM
ID=$AD
IH=$AH
IMN=$AMN
ITER=`expr $ITER + 1`
done	### MAIN LOOP ###
echo "NORMAL END"
