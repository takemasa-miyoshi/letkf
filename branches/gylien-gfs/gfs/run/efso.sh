#!/bin/bash
##===============================================================================
#
#  Run EFSO (ensemble-based forecast sensitivity to observation) for GFS-LETKF
#  Created  August   2013, Daisuke Hotta
#  Modified December 2013, Guo-Yuan Lien
#
#===============================================================================

RUNDIR="$( cd "$( dirname "$0" )" && pwd )"
cd "$RUNDIR"
if [ -f configure.sh ]; then
  . configure.sh
else
  echo "[Error] $0: 'configure.sh' does not exist." 1>&2
  exit 1
fi
. distribute.sh
. stageinout.sh
. datetime.sh

nparts=3
partname[1]='Copy necessary files (stagein)'
partname[2]='Run EFSO'
partname[3]='Collect outputs (stageout)'

if [ "$#" -lt 1 ]; then
  cat 1>&2 << EOF

[efso.sh] Run an EFSO for GFS-LETKF.
           *use settings in 'configure.sh'

Usage: $0 TIME [EFT] [LOCADV_RATE] [WMOIST] [LON1] [LON2] [LAT1] [LAT2] [LEV1] [LEV2]

  TIME         Valid time for EFSO estimation (format: YYYYMMDDHH)
  EFT          Evaluation forecast time (hours)
               (default: 24)
  LOCADV_RATE  Localization advection rate relative to the phase velocity (winds)
               0: No advection
               (default: 0)
  WMOIST       Wight for the moist term in the energy norm (dimensionless)
               (default: 1)
  LON1         Lowest longitude of the target domain (degrees East)
               (default: 0)
  LON2         Higest longitude of the target domain (degrees East)
               (default: 360)
  LAT1         Lowest latitude of the target domain (degrees North)
               (default: -90)
  LAT2         Higest latitude of the target domain (degrees North)
               (default: 90)
  LEV1         Lowest model level of the target domain (level)
               (default: 1)
  LEV2         Higest model level of the target domain (level)
               (default: 64)

EOF
  exit 1
fi

#-------------------------------------------------------------------------------

TIME=$(datetime $1)
yyyymmddhh=${TIME:0:10}
EFT=${2:-24}
LOCADV_RATE="${3:-0}"
WMOIST=${4:-1}
LON1=${5:-0}
LON2=${6:-360}
LAT1=${7:--90}
LAT2=${8:-90}
LEV1=${9:-1}
LEV2=${10:-64}

EVTIME=$(datetime $TIME $EFT h)
EVyyyymmddhh=${EVTIME:0:10}
M6TIME=$(datetime $TIME -6 h)
M6yyyymmddhh=${M6TIME:0:10}

#-------------------------------------------------------------------------------

function now_format {
  date +'%Y-%m-%d %H:%M:%S'
}

function part_timing {
  local p="$1"
  echo "[$(now_format)] ${partname[$p]}" 1>&2
  echo
  printf " %2d. %-55s\n" ${p} "${partname[$p]}"
}

#-------------------------------------------------------------------------------

mkdir -p $TMPMPI
tmpnode="$TMPMPI/node"
tmpstagein="$TMPMPI/stagein"
tmpstageout="$TMPMPI/stageout"

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_efso"
if [ "$SHAREDISK" = '0' ]; then
  ltmprun1="$LTMP1/${tmpsubdir}/run"
  ltmprun2="$LTMP2/${tmpsubdir}/run"
  ltmpout="$LTMP1/${tmpsubdir}/out"
elif [ "$SHAREDISK" = '1' ]; then
  ltmprun1="$OUTDIR/tmp"
  ltmprun2="$OUTDIR/tmp"
  ltmpout="$OUTDIR"
else
  echo "[Error] $0: Unsupported \$SHAREDISK setting." 1>&2
  exit 1
fi

ltmpprog="$ltmprun1/prog"
ltmpobs="$ltmprun1/obs"
ltmpanal="$ltmprun1/anal"
ltmpefso="$ltmprun1/efso"

#ltmpfix="$ltmprun2/fix"

efsoexec=`printf 'efso%03d' $MEMBER`

echo "[$(now_format)] Start efso.sh $@" 1>&2

#===============================================================================
# Print job information

echo
echo " +----------------------------------------------------------------+"
echo " |                       EFSO on GFS-LETKF                        |"
echo " +----------------------------------------------------------------+"
echo
echo "  Validation time:      ${TIME:0:4}-${TIME:4:2}-${TIME:6:2} ${TIME:8:2}:${TIME:10:2}"
echo "  Evaluation FT:        $EFT h"
echo "  Moist term weight:    $WMOIST"
echo "  Target Domain:        (lon: $LON1, $LON2) x"
echo "                        (lat: $LAT1, $LAT2) x"
echo "                        (lev: $LEV1, $LEV2)"
if [ "$LOCADV_RATE" = '0' ]; then
  echo "  Loc advection rate:   Not advected"
else
  echo "  Loc advection rate:   $LOCADV_RATE times phase velocity (winds)"
fi

echo

distribute_da_cycle machinefile 1 1

echo
echo "===================================================================="

#===============================================================================
# 1. Copy necessary files
p=1
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
mkdir -p $tmpstagein
rm -f $tmpstagein/*

# Executable files
cat > $tmpstagein/run.1 << EOF
ck|${DIR}/letkf/${efsoexec}|prog/${efsoexec}
EOF

# Reference analysis at the evaluation time
cat >> $tmpstagein/run.1 << EOF
ne|$ANLGRD/${EVyyyymmddhh}.grd|anal/EVanal.grd
EOF

if [ "$SHAREDISK" = '0' ]; then
  for m in `seq $MEMBER`; do
# Ensemble forecats at the evaluation time
    echo "ck|fcstg/${yyyymmddhh}/${name_m[$m]}/${EVyyyymmddhh}.grd" >> $tmpstagein/out.${node_m[$m]}
# obsanal: observations in ensemble analyses
    echo "ck|obs/obsdiag/${yyyymmddhh}/obsanal_${name_m[$m]}.dat" >> $tmpstagein/out.${node_m[$m]}
  done

# Mean background and analysis
  echo "ck|guesg/mean/${yyyymmddhh}.grd" >> $tmpstagein/out.${node_m[$mmean]}
  echo "ck|analg/mean/${yyyymmddhh}.grd" >> $tmpstagein/out.${node_m[$mmean]}
# Reference analysis at the evaluation time (self analysis)
#  echo "ck|analg/mean/${EVyyyymmddhh}.grd" >> $tmpstagein/out.${node_m[$mmean]}
# Mean forecasts (from 0h and -6h) at the evaluation time
  echo "ck|fcstg/${yyyymmddhh}/mean/${EVyyyymmddhh}.grd" >> $tmpstagein/out.${node_m[$mmean]}
  echo "ck|fcstg/${M6yyyymmddhh}/mean/${EVyyyymmddhh}.grd" >> $tmpstagein/out.${node_m[$mmean]}
# obsgues: observations in mean background
  echo "ck|obs/obsdiag/${yyyymmddhh}/obsgues_mean.dat" >> $tmpstagein/out.${node_m[$mmean]}
fi

stagein $ltmprun1 $ltmprun2 $ltmpout 1  # clean stagein

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
#===============================================================================
# 2. EFSO
p=2
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
if [ "$SHAREDISK" = '0' ]; then
  $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes \
                  bash -c "mkdir -p $ltmpefso; rm -fr $ltmpefso/*" &
  wait
else
  mkdir -p $ltmpefso
  rm -fr $ltmpefso/*
fi

#-------------------------------------------------------------------------------

cat > efso_21.sh << EOF
mem="\$1"
cd $ltmpefso
ln -fs $ltmpout/fcstg/${yyyymmddhh}/\${mem}/${EVyyyymmddhh}.grd fc01\${mem}.grd
ln -fs $ltmpout/obs/obsdiag/${yyyymmddhh}/obsanal_\${mem}.dat obsanal\${mem}.dat
EOF

cat >> efso_21_node1.sh <<EOF
cd $ltmpefso
ln -fs $ltmpout/guesg/mean/${yyyymmddhh}.grd gmean.grd
ln -fs $ltmpout/analg/mean/${yyyymmddhh}.grd anal0.grd
ln -fs $ltmpanal/EVanal.grd amean.grd
#ln -fs $ltmpout/analg/mean/${EVyyyymmddhh}.grd amean.grd
ln -fs $ltmpout/fcstg/${yyyymmddhh}/mean/${EVyyyymmddhh}.grd fme00.grd
ln -fs $ltmpout/fcstg/${M6yyyymmddhh}/mean/${EVyyyymmddhh}.grd fme06.grd
ln -fs $ltmpout/obs/obsdiag/${yyyymmddhh}/obsgues_mean.dat obsgues_me.dat
EOF

#-------------------------------------------------------------------------------

echo
echo "  Prepare EFSO"
echo
ppnl=$((ppn*2+2))
np=1
pcount=0
for m in `seq $MEMBER`; do
  if [ "${node_m[$m]}" = "${node[1]}" ]; then
    pcount=$((pcount+np))
    if [ "$pcount" -gt "$ppnl" ]; then
      echo "    wait..."
      wait
      pcount=$np
    fi
  fi
  echo "  member ${name_m[$m]} on node ${node_m[$m]}"
  $MPIBIN/mpiexec -host ${node_m[$m]} bash efso_21.sh ${name_m[$m]} &
  sleep 0.05s
done

pcount=$((pcount+np))
if [ "$pcount" -gt "$ppnl" ]; then
  echo "    wait..."
  wait
fi
echo "  mean       on node ${node_m[$mmean]}"
$MPIBIN/mpiexec -host ${node_m[$mmean]} bash efso_21_node1.sh &
echo "    wait..."
wait

#-------------------------------------------------------------------------------

cat > efso_22.sh << EOF
cd $ltmpefso
ln -fs $ltmpprog/$efsoexec .

cat > namelist.efso << EOG 
&EFSOPRM
 wmoist=$WMOIST,
 EFT=$EFT,
 locadv_rate=$LOCADV_RATE,
 tar_minlon=$LON1,
 tar_maxlon=$LON2,
 tar_minlat=$LAT1,
 tar_maxlat=$LAT2,
 tar_minlev=$LEV1,
 tar_maxlev=$LEV2,
/
EOG
EOF

if [ "$SHAREDISK" = '0' ]; then
  $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes bash efso_22.sh &
  wait
else
  bash efso_22.sh &
  wait
fi

#-------------------------------------------------------------------------------

np=`cat $tmpnode/machinefile | wc -l`
echo
echo "  Run EFSO on $np processor cores..."

$MPIBIN/mpiexec -machinefile $tmpnode/machinefile -n $np \
                -wdir $ltmpefso ./$efsoexec > /dev/null 2>&1 &
wait
echo "  Done!"

#-------------------------------------------------------------------------------

cat > efso_23_node1.sh << EOF
mkdir -p $ltmpout/efso/ke
mkdir -p $ltmpout/efso/pe
mkdir -p $ltmpout/efso/me
mv -f $ltmpefso/osenseKE.dat $ltmpout/efso/ke/osense_${yyyymmddhh}.dat
mv -f $ltmpefso/osensePE.dat $ltmpout/efso/pe/osense_${yyyymmddhh}.dat
mv -f $ltmpefso/osenseME.dat $ltmpout/efso/me/osense_${yyyymmddhh}.dat
mkdir -p $ltmpout/log/efso
mv -f $ltmpefso/NOUT-000 $ltmpout/log/efso/efso_${yyyymmddhh}.log
EOF

$MPIBIN/mpiexec -host ${node_m[$mmean]} bash efso_23_node1.sh &
wait

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
#===============================================================================
# 3. Collect outputs (stageout)
p=3
part_timing $p
#-------------------------------------------------------------------------------
if [ "$SHAREDISK" = '0' ]; then
#-------------------------------------------------------------------------------

cd $TMPMPI
mkdir -p $tmpstageout
rm -f $tmpstageout/*

for m in `seq $MEMBER`; do
  echo "rm|fcstg/${yyyymmddhh}/${name_m[$m]}/${EVyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$m]}
  echo "rm|obs/obsdiag/${yyyymmddhh}/obsanal_${name_m[$m]}.dat" >> $tmpstageout/out.${node_m[$m]}
done
cat >> $tmpstageout/out.${node_m[$mmean]} << EOF
rm|guesg/mean/${yyyymmddhh}.grd
rm|analg/mean/${yyyymmddhh}.grd
rm|fcstg/${yyyymmddhh}/mean/${EVyyyymmddhh}.grd
rm|fcstg/${M6yyyymmddhh}/mean/${EVyyyymmddhh}.grd
rm|obs/obsdiag/${yyyymmddhh}/obsgues_mean.dat
mv|efso/ke/osense_${yyyymmddhh}.dat
mv|efso/pe/osense_${yyyymmddhh}.dat
mv|efso/me/osense_${yyyymmddhh}.dat
mv|log/efso/efso_${yyyymmddhh}.log
EOF

#-------------------------------------------------------------------------------

stageout $ltmpout 0  # clean stageout
$MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes \
                rm -fr $LTMP1/${tmpsubdir} $LTMP2/${tmpsubdir} &
wait

#-------------------------------------------------------------------------------
elif [ "$SHAREDISK" = '1' ]; then
#-------------------------------------------------------------------------------

rm -fr $ltmprun1 $ltmprun2

#-------------------------------------------------------------------------------
fi
#-------------------------------------------------------------------------------
echo
echo "===================================================================="
#===============================================================================

#rm -fr $TMPMPI

echo
echo " +----------------------------------------------------------------+"
echo " |            EFSO on GFS-LETKF successfully completed            |"
echo " +----------------------------------------------------------------+"
echo

echo "[$(now_format)] Finish efso.sh $@" 1>&2
exit 0

#===============================================================================
