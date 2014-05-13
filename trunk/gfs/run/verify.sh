#!/bin/bash
#===============================================================================
#
#  Run ensemble verifications.
#  Created  November  2012, Guo-Yuan Lien
#  Modified May       2014, Guo-Yuan Lien
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
partname[2]='Run verification'
partname[3]='Collect outputs (stageout)'

if [ "$#" -lt 1 ]; then
  cat 1>&2 << EOF

[verify.sh] Compute verification of ensemble forecasts.
            *use settings in 'configure.sh'

Usage: $0 VTIME [MEMBERS]

  VTIME    Verification time (format: YYYYMMDDHH)
  MEMBERS  List of forecast members ('mean' for ensemble mean)
           'all':   run all members including ensemble mean
           'mems':  run all members but not including ensemble mean
           '2 4 6': run members 2, 4, 6
           (default: 'all')

EOF
  exit 1
fi

#-------------------------------------------------------------------------------

VTIME=$(datetime $1)
Vyyyymmddhh=${VTIME:0:10}
MEMBERS="${2:-all}"

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

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_verify"
if [ "$SHAREDISK" = '0' ]; then
  ltmprun1="$LTMP1/${tmpsubdir}/run"
  ltmprun2="$LTMP2/${tmpsubdir}/run"
  ltmpout="$LTMP1/${tmpsubdir}/out"
  echo "[Warning] $0: It is strongly recommended to use a shared disk (\$SHAREDISK = 1) for verification computation." 1>&2
  echo "[Warning] $0: \$OBS, \$ANLGRDP, \$ANLGRDP2 need to be in a shared disk that can be seen from all computing nodes." 1>&2
elif [ "$SHAREDISK" = '1' ]; then
  ltmprun1="$OUTDIR/tmp"
  ltmprun2="$OUTDIR/tmp"
  ltmpout="$OUTDIR"
else
  echo "[Error] $0: Unsupported \$SHAREDISK setting." 1>&2
  exit 1
fi

ltmpprog="$ltmprun1/prog"
ltmpverify="$ltmprun1/verify"

echo "[$(now_format)] Start verify.sh $@" 1>&2

#-------------------------------------------------------------------------------

if [ "$MEMBERS" = 'all' ] || [ "$MEMBERS" = 'mems' ]; then
  if [ "$MEMBERS" = 'all' ]; then
    MEMBERS='mean '
  elif [ "$MEMBERS" = 'mems' ]; then
    MEMBERS=''
  fi
  for m in `seq $MEMBER`; do
    name_m=`printf '%03d' $m`
    MEMBERS="${MEMBERS}${name_m} "
  done
fi

#===============================================================================
# Print job information

echo
echo " +----------------------------------------------------------------+"
echo " |                          Verification                          |"
echo " +----------------------------------------------------------------+"
echo
echo "  Verification time:    ${VTIME:0:4}-${VTIME:4:2}-${VTIME:6:2} ${VTIME:8:2}:${VTIME:10:2}"
echo "  Forecast length:      $FCSTLEN h"
echo "  Output interval:      $FCSTOUT h"
echo

distribute_fcst machinefile "$MEMBERS" 1 1 1

echo
echo "  # verification members: $fmember"
echo
for m in `seq $fmember`; do
  if [ "$SHAREDISK" = '0' ]; then
    echo "    ${name_m[$m]} stored on ${node_m[$m]}"
  else
    echo "    ${name_m[$m]} run on ${nodes_gfs_m[$m]}"
  fi
done
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
ck|${DIR}/verify/verify|prog/verify
EOF

# Ensemble forecast files
if [ "$SHAREDISK" = '0' ]; then
  for m in `seq $fmember`; do
    fh=0
    while [ "$fh" -le "$FCSTLEN" ]; do
      Syyyymmddhh=$(datetime $VTIME -$fh h | cut -c 1-10)
      if [ -s "$OUTDIR/fcstg/${Syyyymmddhh}/${name_m[$m]}/${Vyyyymmddhh}.grd" ] &&
         [ -s "$OUTDIR/fcstgp/${Syyyymmddhh}/${name_m[$m]}/${Vyyyymmddhh}.grd" ]; then
        echo "ne|fcstg/${Syyyymmddhh}/${name_m[$m]}/${Vyyyymmddhh}.grd" >> $tmpstagein/out.${node_m[$m]}
        echo "ne|fcstgp/${Syyyymmddhh}/${name_m[$m]}/${Vyyyymmddhh}.grd" >> $tmpstagein/out.${node_m[$m]}
      fi
    fh=$((fh+FCSTOUT))
    done
  done
fi

stagein $ltmprun1 $ltmprun2 $ltmpout 1  # clean stagein

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
#===============================================================================
# 2. Verification
p=2
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
cp -f $RUNDIR/datetime.sh .

cat > verify_21.sh << EOF
. datetime.sh
mem="\$1"

mkdir -p $ltmpverify/\${mem}
cd $ltmpverify/\${mem}
ln -fs $ltmpprog/verify .
fh=0
while [ "\$fh" -le "$FCSTLEN" ]; do
  fhhh=\$(printf '%03d' \$fh)
  Syyyymmddhh=\$(datetime $VTIME -\$fh h | cut -c 1-10)

  if [ -s "$ltmpout/fcstg/\${Syyyymmddhh}/\${mem}/${Vyyyymmddhh}.grd" ] &&
     [ -s "$ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}/${Vyyyymmddhh}.grd" ]; then
    cd $ltmpverify/\${mem}
    rm -f fcst.grd fcstp.grd obs??.dat ana??.grd
    ln -s $ltmpout/fcstg/\${Syyyymmddhh}/\${mem}/${Vyyyymmddhh}.grd fcst.grd
    ln -s $ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}/${Vyyyymmddhh}.grd fcstp.grd
###### only support shared disk
    cat $OBS/obs${Vyyyymmddhh}/t.dat > obs01.dat
    cat $OBS/obs${Vyyyymmddhh}/t-1.dat >> obs01.dat
    cat $OBS/obs${Vyyyymmddhh}/t+1.dat >> obs01.dat
    ln -s $ANLGRDP/${Vyyyymmddhh}.grd ana01.grd
#    ln -s $ANLGRDP2/${Vyyyymmddhh}.grd ana02.grd
######
    ./verify > /dev/null 2>&1

    mkdir -p $ltmpout/verfo1/\${fhhh}/\${mem}
    mkdir -p $ltmpout/verfa1/\${fhhh}/\${mem}
    mkdir -p $ltmpout/verfa2/\${fhhh}/\${mem}
    mv -f vrfobs01.dat $ltmpout/verfo1/\${fhhh}/\${mem}/${Vyyyymmddhh}.dat
    mv -f vrfana01.dat $ltmpout/verfa1/\${fhhh}/\${mem}/${Vyyyymmddhh}.dat
#    mv -f vrfana02.dat $ltmpout/verfa2/\${fhhh}/\${mem}/${Vyyyymmddhh}.dat
  fi
fh=\$((fh+$FCSTOUT))
done
EOF

#-------------------------------------------------------------------------------

echo
ppnl=$ppn
#ppnl=$((ppn*2))
np=1
pcount=0
for m in `seq $fmember`; do
  if [ "${node_m[$m]}" = "${node[1]}" ]; then
    pcount=$((pcount+np))
    if [ "$pcount" -gt "$ppnl" ]; then
      echo "    wait..."
      wait
      pcount=$np
    fi
  fi
  echo "  member ${name_m[$m]} on node ${node_m[$m]}"
  $MPIBIN/mpiexec -host ${node_m[$m]} bash verify_21.sh ${name_m[$m]} &
  sleep 0.05s
done
echo "    wait..."
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

for m in `seq $fmember`; do
  fh=0
  while [ "$fh" -le "$FCSTLEN" ]; do
    fhhh=`printf '%03d' $fh`
    Syyyymmddhh=$(datetime $VTIME -$fh h | cut -c 1-10)
    if [ -s "$OUTDIR/fcstg/${Syyyymmddhh}/${name_m[$m]}/${Vyyyymmddhh}.grd" ] &&
       [ -s "$OUTDIR/fcstgp/${Syyyymmddhh}/${name_m[$m]}/${Vyyyymmddhh}.grd" ]; then
      echo "rm|fcstg/${Syyyymmddhh}/${name_m[$m]}/${Vyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$m]}
      echo "rm|fcstgp/${Syyyymmddhh}/${name_m[$m]}/${Vyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$m]}
      echo "mv|verfo1/${fhhh}/${name_m[$m]}/${Vyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$m]}
      echo "mv|verfa1/${fhhh}/${name_m[$m]}/${Vyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$m]}
      echo "mv|verfa2/${fhhh}/${name_m[$m]}/${Vyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$m]}
    fi
  fh=$((fh+FCSTOUT))
  done
done

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
echo " |              Verification successfully completed               |"
echo " +----------------------------------------------------------------+"
echo

echo "[$(now_format)] Finish verify.sh $@" 1>&2
exit 0

#===============================================================================
