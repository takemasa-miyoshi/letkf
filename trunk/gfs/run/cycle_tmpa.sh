#!/bin/bash
#===============================================================================
#
#  Run one GFS-LETKF cycle.
#  Created  October   2012, Guo-Yuan Lien
#  Modified April     2013, Guo-Yuan Lien
#  Modified September 2013, Guo-Yuan Lien
#  Modified December  2013, Guo-Yuan Lien
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

nparts=10
partname[1]='Copy necessary files (stagein)'
partname[2]='Run ensemble forecasts'
partname[3]='Compute ensemble mean'
partname[4]='Convert GFS sig/sfc format into grid format for LETKF'
if [ "$THIN_OPT" = '2' ]; then
  partname[5]='Run superobing/thinning'
elif [ "$OBSOPE_OPT" = '2' ]; then
  partname[5]='Run obs operator for ens mean (using GSI, may do thinning)'
else
  partname[5]=' -- skipped with current settings'
fi
if [ "$OBSOPE_OPT" = '1' ]; then
  partname[6]=' -- skipped with current settings'
elif [ "$OBSOPE_OPT" = '2' ]; then
  partname[6]='Run obs operator for each member (using GSI)'
else
  echo "[Error] $0: Unsupported \$OBSOPE_OPT setting." 1>&2
  exit 1
fi
if [ "$OBSOPE_OPT" = '1' ]; then
  partname[7]='Run obs operator for each member (using obsope)'
elif [ "$OBSOPE_OPT" = '2' ]; then
  partname[7]='Convert GSI diag files into LETKF obs2 format'
fi
partname[8]='Run LETKF'
partname[9]='Convert LETKF output into GFS sig/sfc format'
partname[10]='Collect outputs (stageout)'

#-------------------------------------------------------------------------------

if [ "$#" -lt 1 ]; then
  cat 1>&2 << EOF

[cycle.sh] Run a forecast/analysis cycle for GFS-LETKF.
           *use settings in 'configure.sh'

Usage: $0 STIME [PART] [MCYCLE]

  STIME   Start time of the cycle (format: YYYYMMDDHH)
  PART    Parts of the script to be run
EOF
for p in `seq $nparts`; do
  echo "          ${p}:     ${partname[$p]}" 1>&2
done
  cat 1>&2 << EOF
          'all': Run all parts
          '3-6': Run parts 3 to 6
          '4-' : Run parts after 4
          (default: all)
  MCYCLE  Options for multiple cycle run
          0: Not a multiple cycle run (clean stagein/out)
          1: The first cycle of a multiple cycle run (clean stagein)
          2: A middle cycle of a multiple cycle run
          3: The Last cycle of a multiple cycle run (clean stageout)
          (not clean stagein  - trust local files from the previous cycle)
          (not clean stageout - keep local files for the next cycle)
          (default: 0)

EOF
  exit 1
fi

#-------------------------------------------------------------------------------

STIME=$(datetime $1)
ETIME=$(datetime $STIME $LCYCLE h)
Syyyymmddhh=${STIME:0:10}
Eyyyymmddhh=${ETIME:0:10}
EEyyyymmddhh=$(datetime $ETIME $LCYCLE h | cut -c 1-10)
PART=${2:-all}
MCYCLE=${3:-0}

if [ "$ADAPTINFL" -eq 1 ]; then
  INFLyyyymmddhh=${Eyyyymmddhh}
elif [ "$ADAPTINFL" -eq 2 ]; then
  INFLyyyymmddhh=${Syyyymmddhh}
fi

#-------------------------------------------------------------------------------

if [ "$PART" = 'all' ]; then
  p1=1
  p2=$nparts
elif [ `echo $PART | grep '-'` ]; then
  p1=`echo $PART | cut -d '-' -s -f1`
  if [ -z "$p1" ]; then
    p1=1
  fi
  p2=`echo $PART | cut -d '-' -s -f2`
  if [ -z "$p2" ]; then
    p2=$nparts
  fi
else
  p1=$PART
  p2=$PART
fi
for p in `seq $nparts`; do
  if [ "$p" -ge "$p1" ] && [ "$p" -le "$p2" ]; then
    run_part[$p]=1
  else
    run_part[$p]=0
  fi
done

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
tmpdata="$TMPMPI/data"
tmpstagein="$TMPMPI/stagein"
tmpstageout="$TMPMPI/stageout"

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_cycle"
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
ltmpssio="$ltmprun1/ssio"
ltmpgfs="$ltmprun1/gfs"
ltmpsuperob="$ltmprun1/superob"
ltmpobsope="$ltmprun1/obsope"
ltmpgsi="$ltmprun1/gsi"
ltmpreaddiag="$ltmprun1/readdiag"
ltmpletkf="$ltmprun1/letkf"

ltmpfix="$ltmprun2/fix"
ltmpfixgsi="$ltmprun2/fixgsi"
ltmpfixcrtm="$ltmprun2/fixcrtm"

letkfexec=`printf 'letkf%03d' $MEMBER`
gfsmeanexec=`printf 'gfsmean%03d' $MEMBER`

#-------------------------------------------------------------------------------

if [ ! -f "$DIR/letkf/${letkfexec}" ]; then
  echo "[Error] $0: '$DIR/letkf/${letkfexec}' does not exist." 1>&2
  exit 1
fi
if [ ! -f "$DIR/util/${gfsmeanexec}" ]; then
  echo "[Error] $0: '$DIR/util/${gfsmeanexec}' does not exist." 1>&2
  exit 1
fi
if [ "$OBSOPE_OPT" = '2' ]; then
  gsiwindow=$((WINDOW_E-LCYCLE))
  gsiwindow2=$((LCYCLE-WINDOW_S))
  if [ "$gsiwindow" -ne "$gsiwindow2" ]; then
    echo "[Error] $0: 'Assimilation window needs to be symmetric when using GSI." 1>&2
    exit 1
  fi
fi

echo "[$(now_format)] Start cycye.sh $@" 1>&2

#-------------------------------------------------------------------------------

gsi_obstypes='conv'
#gsi_obstypes='hirs2_n14 msu_n14 sndr_g08 sndr_g11 sndr_g11 sndr_g12 sndr_g13 sndr_g08_prep sndr_g11_prep sndr_g12_prep sndr_g13_prep sndrd1_g11 sndrd2_g11 sndrd3_g11 sndrd4_g11 sndrd1_g12 sndrd2_g12 sndrd3_g12 sndrd4_g12 sndrd1_g13 sndrd2_g13 sndrd3_g13 sndrd4_g13 hirs3_n15 hirs3_n16 hirs3_n17 amsua_n15 amsua_n16 amsua_n17 amsub_n15 amsub_n16 amsub_n17 hsb_aqua airs_aqua amsua_aqua imgr_g08 imgr_g11 imgr_g12 pcp_ssmi_dmsp pcp_tmi_trmm conv sbuv2_n16 sbuv2_n17 sbuv2_n18 sbuv2_n19 gome_metop-a omi_aura ssmi_f13 ssmi_f14 ssmi_f15 hirs4_n18 hirs4_metop-a amsua_n18 amsua_metop-a mhs_n18 mhs_metop-a amsre_low_aqua amsre_mid_aqua amsre_hig_aqua ssmis_f16 ssmis_f17 ssmis_f18 iasi_metop-a hirs4_n19 amsua_n19 mhs_n19 seviri_m08 seviri_m09 seviri_m10 atms_npp cris_npp'

#-------------------------------------------------------------------------------

time=$(datetime $STIME $WINDOW_S h)
is=0
while [ $(datetime $time) -le $(datetime $STIME $WINDOW_E h) ]; do
  is=$((is+1))
  time_sl[$is]="${time:0:4}-${time:4:2}-${time:6:2} ${time:8:2}:${time:10:2}"
  fh=$((WINDOW_S+LTIMESLOT*(is-1)))
  fhh_sl[$is]=`printf '%02d' $fh`
  if [ "$fh" -eq "$LCYCLE" ]; then
    baseslot=$is
    obsfile_sl[$is]="t.dat"
  elif [ "$fh" -lt "$LCYCLE" ]; then
    obsfile_sl[$is]="t-$((LCYCLE-fh)).dat"
  else
    obsfile_sl[$is]="t+$((fh-LCYCLE)).dat"
  fi
time=$(datetime $time $LTIMESLOT h)
done
nslots=$is

time=$(datetime $STIME $WINDOW_S h)
is=0
while [ $(datetime $time) -le $(datetime $STIME $WINDOW_E h) ]; do
  is=$((is+1))
  time_slmean[$is]="${time:0:4}-${time:4:2}-${time:6:2} ${time:8:2}:${time:10:2}"
  fh=$((WINDOW_S+LTIMESLOTMEAN*(is-1)))
  fhh_slmean[$is]=`printf '%02d' $fh`
time=$(datetime $time $LTIMESLOTMEAN h)
done
nslotsmean=$is

#===============================================================================
# Print job information

echo
echo " +----------------------------------------------------------------+"
echo " |                           GFS-LETKF                            |"
echo " +----------------------------------------------------------------+"
for p in `seq $nparts`; do
  if [ "${run_part[$p]}" = '1' ]; then
    printf " | %2d. %-58s |\n" ${p} "${partname[$p]}"
  fi
done
echo " +----------------------------------------------------------------+"
echo
echo "  State time:           ${STIME:0:4}-${STIME:4:2}-${STIME:6:2} ${STIME:8:2}:${STIME:10:2}"
echo "  Forecast length:      $WINDOW_E h"
echo "  Assimilation window:  $((WINDOW_E-WINDOW_S)) h"
echo
echo "  Observation timeslots:"
is=1
while [ "$is" -le "$nslots" ]; do
  isf=`printf '%02d' $is`
  echo "    $isf - ${time_sl[$is]}"
is=$((is+1))
done
echo

distribute_da_cycle machinefile 1 1

echo
echo "===================================================================="

#===============================================================================
# 1. Copy necessary files
p=1
if [ "${run_part[$p]}" = '1' ]; then
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
mkdir -p $tmpstagein
rm -f $tmpstagein/*

# Fix files
cat > $tmpstagein/run.2 << EOF
ck|${FIXGLOBAL}|fix
ck|${FIXGSI}|fixgsi
ck|${FIXCRTM}|fixcrtm
EOF

# Executable files
cat > $tmpstagein/run.1 << EOF
ck|${EXECGLOBAL}/global_fcst|prog/global_fcst
ck|${EXECGLOBAL}/global_gsi|prog/global_gsi
ck|${EXECGLOBAL}/global_sfchdr|prog/global_sfchdr
ck|${EXECGLOBAL}/global_sighdr|prog/global_sighdr
ck|${DIR}/letkf/${letkfexec}|prog/${letkfexec}
ck|${DIR}/letkf/obsope|prog/obsope
ck|${DIR}/obs/superob|prog/superob
ck|${DIR}/obs/readdiag_conv|prog/readdiag_conv
ck|${DIR}/ssio/grd2ss|prog/grd2ss
ck|${DIR}/ssio/ss2grd|prog/ss2grd
ck|${DIR}/ssio/ss2grdp|prog/ss2grdp
ck|${DIR}/util/${gfsmeanexec}|prog/${gfsmeanexec}
ck|${DIR}/verify/verify|prog/verify
EOF

# Reference analysis
cat >> $tmpstagein/run.1 << EOF
ne|$ANLGFS/${Eyyyymmddhh}.sig|anal/Eanal.sig
ne|$ANLGFS/${Eyyyymmddhh}.sfc|anal/Eanal.sfc
EOF

# Observations
if [ "$OBSOPE_OPT" = '1' ]; then
  for is in `seq $nslots`; do
    isf=`printf '%02d' $is`
    echo "ne|$OBS/obs${Eyyyymmddhh}/${obsfile_sl[$is]}|obs/obs${isf}.dat" \
         >> $tmpstagein/run.1
  done
  echo "ne|$OBSPP/tmpa_obs${Eyyyymmddhh}.dat|obs/obsmn.dat" >> $tmpstagein/run.1

  Emm=$(datetime $ETIME | cut -c 5-6)
  Edd=$(datetime $ETIME | cut -c 7-8)
  nEmm=$((10#${Emm}))
  nEdd=$((10#${Edd}))
  if [ "$nEdd" -le 10 ]; then
    tpm=1
  elif [ "$nEdd" -le 20 ]; then
    tpm=2
  else
    tpm=3
  fi
  tp=$(((nEmm-1)*3+tpm))
  tpf=`printf '%03d' $tp`
  echo "ne|$PPCDFM/cdf-${tpf}.dat|obs/cdfm.grd" >> $tmpstagein/run.1
  echo "ne|$PPCDFO/cdf-${tpf}.dat|obs/cdfo.grd" >> $tmpstagein/run.1
  if [ -s "$PPCORR/corr-${tpf}.dat" ]; then
    echo "ne|$PPCORR/corr-${tpf}.dat|obs/ppmask.grd" >> $tmpstagein/run.1
  fi
elif [ "$OBSOPE_OPT" = '2' ]; then
  cat >> $tmpstagein/run.1 << EOF
ne|$OBSNCEP/obs${Eyyyymmddhh}/prepbufr.gdas.${Eyyyymmddhh}.nr|obs/prepbufr.gdas.nr
EOF
fi

if [ "$SHAREDISK" = '0' ]; then
# Ensemble analyses from the previous cycle
  for m in `seq $MEMBER`; do
    echo "ck|anal/${name_m[$m]}/${Syyyymmddhh}.sig" >> $tmpstagein/out.${node_m[$m]}
    echo "ck|anal/${name_m[$m]}/${Syyyymmddhh}.sfc" >> $tmpstagein/out.${node_m[$m]}
  done

# Adaptive inflation field from the previous cycle
  if [ "$ADAPTINFL" -ge 1 ] && [ -s "$OUTDIR/infl/${INFLyyyymmddhh}.grd" ]; then
    echo "ck|infl/${INFLyyyymmddhh}.grd" >> $tmpstagein/out.${node_m[$mmean]}
  fi
fi

#-------------------------------------------------------------------------------

if [ "$MCYCLE" = '0' ] || [ "$MCYCLE" = '1' ]; then
  stagein $ltmprun1 $ltmprun2 $ltmpout 1  # clean stagein
else
  stagein $ltmprun1 $ltmprun2 $ltmpout 0  # not clean stagein, trust local files from the previous cycle
fi

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
#===============================================================================
# 2. Ensemble forecasts
p=2
if [ "${run_part[$p]}" = '1' ]; then
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
cp -f $RUNDIR/run_gfs.sh .
for m in `seq $MEMBER`; do
  $MPIBIN/mpiexec -host ${node_m[$m]} bash run_gfs.sh \
                  $ltmpout/anal/${name_m[$m]}/${Syyyymmddhh}.sig \
                  $ltmpout/anal/${name_m[$m]}/${Syyyymmddhh}.sfc \
                  $FHMAX $FHOUT $ltmpgfs/${name_m[$m]} \
                  $ltmpprog $ltmpfix &
  sleep 0.05s
done
wait

#-------------------------------------------------------------------------------

if [ "$LOG_OPT" -le 2 ]; then
  mkdir -p $OUTDIR/log/gfs/${Syyyymmddhh}
fi
echo
ppnl=$((ppn*mem_nodes_gfs))
pcount=0
for m in `seq $MEMBER`; do
  np=`cat $tmpnode/machinefile.gfs.${name_m[$m]} | wc -l`
  if [ "${node_m[$m]}" = "${node[1]}" ]; then
    pcount=$((pcount+np))
    if [ "$pcount" -gt "$ppnl" ]; then
      echo "    wait..."
      wait
      pcount=$np
    fi
  fi
  echo "  run GFS for member ${name_m[$m]} on node ${nodes_gfs_m[$m]}"
  if [ "$LOG_OPT" -le 2 ]; then
    $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.gfs.${name_m[$m]} -n $np \
                    -wdir $ltmpgfs/${name_m[$m]} \
                    ./global_fcst > $OUTDIR/log/gfs/${Syyyymmddhh}/gfs_${name_m[$m]}.log 2>&1 &
  else
    $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.gfs.${name_m[$m]} -n $np \
                    -wdir $ltmpgfs/${name_m[$m]} \
                    ./global_fcst > /dev/null 2>&1 &
  fi
  sleep 0.05s
done
echo "    wait..."
wait

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
#===============================================================================
# 3. Compute ensemble mean
p=3
if [ "${run_part[$p]}" = '1' ]; then
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI

if [ "$SHAREDISK" = '0' ]; then
  cat > cycle_31.sh << EOF
fhh="\$1"
EOF
  for m in `seq $MEMBER`; do
    echo "node_m[$m]='${node_m[$m]}'" >> cycle_31.sh
  done
  cat >> cycle_31.sh << EOF
mkdir -p $ltmpgfs/mean
cd $ltmpgfs/mean
rm -f sig* sfc*
ln -fs $ltmpprog/$gfsmeanexec .
prefixlen="\$(expr length \$(hostname))"
for m in \$(seq $MEMBER); do
  mem="\$(printf '%03d' \$m)"
  node_m_tmp="\${node_m[\$m]}"
  if [ "\${node_m_tmp:0:\$prefixlen}" = "\$(hostname)" ]; then
    ln -s $ltmpgfs/\${mem}/SIG.F\${fhh} sig\${mem}
    ln -s $ltmpgfs/\${mem}/SFC.F\${fhh} sfc\${mem}
  fi
done
EOF
else
  mkdir -p $ltmpgfs/mean
  cd $ltmpgfs/mean
  ln -fs $ltmpprog/$gfsmeanexec .
  cd $TMPMPI
fi

cat $tmpnode/machinefile | head -n $MEMBER > cycle_31.machinefile
np=`cat cycle_31.machinefile | wc -l`

#-------------------------------------------------------------------------------
if [ "$OBSOPE_OPT" = '2' ]; then # Compute ens mean at several timeslots (used by GSI)
#-------------------------------------------------------------------------------

  echo
  for is in `seq $nslotsmean`; do
    echo "  compute ensemble mean at ${time_slmean[$is]}"
    if [ "$SHAREDISK" = '0' ]; then
      $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes \
                      bash cycle_31.sh ${fhh_slmean[$is]} &
      wait
    else
      cd $ltmpgfs/mean
      rm -f sig* sfc*
      for m in `seq $MEMBER`; do
        ln -s $ltmpgfs/${name_m[$m]}/SIG.F${fhh_slmean[$is]} sig${name_m[$m]}
        ln -s $ltmpgfs/${name_m[$m]}/SFC.F${fhh_slmean[$is]} sfc${name_m[$m]}
      done
      cd $TMPMPI
    fi

    $MPIBIN/mpiexec -machinefile cycle_31.machinefile -n $np \
                    -wdir $ltmpgfs/mean ./$gfsmeanexec > /dev/null 2>&1 &
    wait
    $MPIBIN/mpiexec -host ${node_m[$mmean]} \
      bash -c "cd $ltmpgfs/mean; mv -f sig_me SIG.F${fhh_slmean[$is]}; mv -f sfc_me SFC.F${fhh_slmean[$is]}" &
    wait
  done

#-------------------------------------------------------------------------------
else # Compute ens mean only at base timeslot
#-------------------------------------------------------------------------------

  echo
  echo "  compute ensemble mean at ${time_sl[$baseslot]}"
  if [ "$SHAREDISK" = '0' ]; then
    $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes \
                    bash cycle_31.sh ${fhh_sl[$baseslot]} &
    wait
  else
    cd $ltmpgfs/mean
    rm -f sig* sfc*
    for m in `seq $MEMBER`; do
      ln -s $ltmpgfs/${name_m[$m]}/SIG.F${fhh_sl[$baseslot]} sig${name_m[$m]}
      ln -s $ltmpgfs/${name_m[$m]}/SFC.F${fhh_sl[$baseslot]} sfc${name_m[$m]}
    done
    cd $TMPMPI
  fi

  $MPIBIN/mpiexec -machinefile cycle_31.machinefile -n $np \
                  -wdir $ltmpgfs/mean ./$gfsmeanexec > /dev/null 2>&1 &
  wait
  $MPIBIN/mpiexec -host ${node_m[$mmean]} \
    bash -c "cd $ltmpgfs/mean; mv -f sig_me SIG.F${fhh_sl[$baseslot]}; mv -f sfc_me SFC.F${fhh_sl[$baseslot]}" &
  wait

#-------------------------------------------------------------------------------
fi
#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
#===============================================================================
# 4. Convert GFS sig/sfc format into grid format for LETKF
p=4
if [ "${run_part[$p]}" = '1' ]; then
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
if [ "$SHAREDISK" = '0' ]; then
  $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes \
                  bash -c "mkdir -p $ltmpletkf; rm -fr $ltmpletkf/*" &
  wait
else
  mkdir -p $ltmpletkf
  rm -fr $ltmpletkf/*
fi

#-------------------------------------------------------------------------------
if [ "$OBSOPE_OPT" = '1' ]; then
#-------------------------------------------------------------------------------

  cat > cycle_41.sh << EOF
mem="\$1"
EOF
for is in `seq $nslots`; do
  echo "fhh_sl[$is]='${fhh_sl[$is]}'" >> cycle_41.sh
done
cat >> cycle_41.sh << EOF
mkdir -p $ltmpobsope/\$mem
rm -rf $ltmpobsope/\${mem}/*
mkdir -p $ltmpssio/\$mem
cd $ltmpssio/\$mem
rm -f fort.*
ln -fs $ltmpprog/ss2grd .
for is in \$(seq $nslots); do
  isf="\$(printf '%02d' \$is)"
  ln -fs $ltmpgfs/\${mem}/SIG.F\${fhh_sl[\$is]} fort.11
  ln -fs $ltmpgfs/\${mem}/SFC.F\${fhh_sl[\$is]} fort.12
  ./ss2grd
  if [ "\$is" -eq "$baseslot" ]; then
    cp -f fort.31 $ltmpletkf/gues\${mem}.grd
EOF
  if [ "$OUT_OPT" -le 2 ]; then
    cat >> cycle_41.sh << EOF
    mkdir -p $ltmpout/guesg/\${mem}
    cp -f fort.31 $ltmpout/guesg/\${mem}/$Eyyyymmddhh.grd
EOF
  fi
  if [ "$OUT_OPT" -le 3 ]; then
    cat >> cycle_41.sh << EOF
    mkdir -p $ltmpout/gues/\${mem}
    cp -f $ltmpgfs/\${mem}/SIG.F\${fhh_sl[\$is]} $ltmpout/gues/\${mem}/$Eyyyymmddhh.sig
    cp -f $ltmpgfs/\${mem}/SFC.F\${fhh_sl[\$is]} $ltmpout/gues/\${mem}/$Eyyyymmddhh.sfc
EOF
  fi
  cat >> cycle_41.sh << EOF
  fi
  mv -f fort.31 $ltmpobsope/\${mem}/gues\${isf}.grd
done
EOF

#-------------------------------------------------------------------------------
elif [ "$OBSOPE_OPT" = '2' ]; then
#-------------------------------------------------------------------------------

  cat > cycle_41.sh << EOF
mem="\$1"
mkdir -p $ltmpssio/\$mem
cd $ltmpssio/\$mem
rm -f fort.*
ln -fs $ltmpprog/ss2grd .
ln -fs $ltmpgfs/\${mem}/SIG.F${fhh_sl[$baseslot]} fort.11
ln -fs $ltmpgfs/\${mem}/SFC.F${fhh_sl[$baseslot]} fort.12
./ss2grd
EOF
  if [ "$OUT_OPT" -le 2 ]; then
    cat >> cycle_41.sh << EOF
mkdir -p $ltmpout/guesg/\${mem}
cp -f fort.31 $ltmpout/guesg/\${mem}/$Eyyyymmddhh.grd
EOF
  fi
  if [ "$OUT_OPT" -le 3 ]; then
    cat >> cycle_41.sh << EOF
mkdir -p $ltmpout/gues/\${mem}
cp -f $ltmpgfs/\${mem}/SIG.F${fhh_sl[$baseslot]} $ltmpout/gues/\${mem}/$Eyyyymmddhh.sig
cp -f $ltmpgfs/\${mem}/SFC.F${fhh_sl[$baseslot]} $ltmpout/gues/\${mem}/$Eyyyymmddhh.sfc
EOF
  fi
  cat >> cycle_41.sh << EOF
mv -f fort.31 $ltmpletkf/gues\${mem}.grd
EOF

#-------------------------------------------------------------------------------
fi
#-------------------------------------------------------------------------------

rm -f cycle_41_node1.sh
if [ "$THIN_OPT" = '2' ]; then
  cat >> cycle_41_node1.sh << EOF
mkdir -p $ltmpssio/mean
cd $ltmpssio/mean
rm -f fort.*
ln -fs $ltmpprog/ss2grd .
ln -fs $ltmpgfs/mean/SIG.F${fhh_sl[$baseslot]} fort.11
ln -fs $ltmpgfs/mean/SFC.F${fhh_sl[$baseslot]} fort.12
./ss2grd
mkdir -p $ltmpsuperob
rm -fr $ltmpsuperob/*
mv -f fort.31 $ltmpsuperob/refpres.dat
EOF
fi
if [ "$OUT_OPT" -le 5 ]; then
  cat >> cycle_41_node1.sh << EOF
mkdir -p $ltmpout/gues/mean
cp -f $ltmpgfs/mean/SIG.F${fhh_sl[$baseslot]} $ltmpout/gues/mean/${Eyyyymmddhh}.sig
cp -f $ltmpgfs/mean/SFC.F${fhh_sl[$baseslot]} $ltmpout/gues/mean/${Eyyyymmddhh}.sfc
EOF
fi

#-------------------------------------------------------------------------------

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
  $MPIBIN/mpiexec -host ${node_m[$m]} bash cycle_41.sh ${name_m[$m]} &
  sleep 0.05s
done

pcount=$((pcount+np))
if [ "$pcount" -gt "$ppnl" ]; then
  echo "    wait..."
  wait
fi
if [ -s "cycle_41_node1.sh" ]; then
  echo "  mean       on node ${node_m[$mmean]}"
  $MPIBIN/mpiexec -host ${node_m[$mmean]} bash cycle_41_node1.sh &
fi
echo "    wait..."
wait

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
#===============================================================================
# 5. (When use obsope) Run superobing/thinning
#    (When use GSI   ) Run obs operator for ensemble mean
p=5
if [ "${run_part[$p]}" = '1' ]; then
if [ "$THIN_OPT" = '2' ] || [ "$OBSOPE_OPT" = '2' ]; then
part_timing $p
#-------------------------------------------------------------------------------
if [ "$THIN_OPT" = '2' ]; then
#-------------------------------------------------------------------------------

cd $TMPMPI
rm -f cycle_51.sh
for is in `seq $nslots`; do
  echo "obsfile_sl[$is]='${obsfile_sl[$is]}'" >> cycle_51.sh
done
cat >> cycle_51.sh << EOF
cd $ltmpsuperob
ln -fs $ltmpprog/superob .
for is in \$(seq $nslots); do
  isf="\$(printf '%02d' \$is)"
  ln -fs $ltmpobs/obs\${isf}.dat .
done
./superob > superob.log 2>&1
for is in \$(seq $nslots); do
  isf="\$(printf '%02d' \$is)"
  rm -f $ltmpobs/obs\${isf}.dat
EOF
if [ "$OBSOUT_OPT" -le 3 ]; then
  cat >> cycle_51.sh << EOF
  mkdir -p $ltmpout/obs/superob/obs${Eyyyymmddhh}
  cp -f sup\${isf}.dat $ltmpout/obs/superob/obs${Eyyyymmddhh}/\${obsfile_sl[\$is]}
EOF
fi
cat >> cycle_51.sh << EOF
  mv -f sup\${isf}.dat $ltmpobs/obs\${isf}.dat
done
EOF
if [ "$LOG_OPT" -le 3 ]; then
  cat >> cycle_51.sh << EOF
mkdir -p $ltmpout/log/superob
mv -f superob.log $ltmpout/log/superob/superob_${Eyyyymmddhh}.log
EOF
fi
if [ "$LOG_OPT" -le 1 ]; then
  cat >> cycle_51.sh << EOF
mkdir -p $ltmpout/log/superob_d
tar czf $ltmpout/log/superob_d/superob_${Eyyyymmddhh}.tar.gz log.*
EOF
fi

echo
echo "  Run superobing/thinning on node ${node_m[$mmean]}..."
$MPIBIN/mpiexec -host ${node_m[$mmean]} bash cycle_51.sh &
wait
echo "  Done!"

#-------------------------------------------------------------------------------
elif [ "$OBSOPE_OPT" = '2' ]; then
#-------------------------------------------------------------------------------

cd $TMPMPI
cp -f $RUNDIR/run_gsi.sh .
$MPIBIN/mpiexec -host ${node_m[$mmean]} bash run_gsi.sh \
                $ltmpgsi/${name_m[$mmean]} $ltmpprog $ltmpfixgsi $ltmpfixcrtm \
                $ltmpgfs/${name_m[$mmean]} $ltmpobs 1 $gsiwindow $THIN_OPT &
wait

echo
echo "  Run GSI for ensemble mean on node ${nodes_gsi_m[$mmean]}..."
np=`cat $tmpnode/machinefile.gsi.${name_m[$mmean]} | wc -l`
if [ "$LOG_OPT" -le 3 ]; then
  mkdir -p $OUTDIR/log/gsi/${Syyyymmddhh}
  $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.gsi.${name_m[$mmean]} -n $np \
                  -wdir $ltmpgsi/${name_m[$mmean]} \
                  ./global_gsi > $OUTDIR/log/gsi/${Syyyymmddhh}/gsi_${name_m[$mmean]}.log 2>&1 &
else
  $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.gsi.${name_m[$mmean]} -n $np \
                  -wdir $ltmpgsi/${name_m[$mmean]} \
                  ./global_gsi > /dev/null 2>&1 &
fi
wait
echo "  Done!"

#-------------------------------------------------------------------------------

cat > cycle_52.sh << EOF
cd $ltmpgsi/${name_m[$mmean]}
EOF
if [ "$OBSOUT_OPT" -le 1 ]; then
  cat >> cycle_52.sh << EOF
mkdir -p $ltmpout/obs/gsidiag/${Eyyyymmddhh}
listall="$gsi_obstypes"
for type in \$listall; do
  count=\$(ls pe*.\${type}_01* 2> /dev/null | wc -l)
  if [ "\$count" -gt 0 ]; then
    cat pe*.\${type}_01* > diag_\${type}
    gzip -c diag_\${type} > $ltmpout/obs/gsidiag/${Eyyyymmddhh}/diag_\${type}_${name_m[$mmean]}.gz
  fi
done
EOF
fi
if [ "$OBSOUT_OPT" -le 3 ]; then
  cat >> cycle_52.sh << EOF
mkdir -p $ltmpout/obs/obsinput/${Eyyyymmddhh}
tar czf $ltmpout/obs/obsinput/${Eyyyymmddhh}/obs_input.tar.gz obs_input.*
EOF
fi
if [ "$LOG_OPT" -le 1 ]; then
  cat >> cycle_52.sh << EOF
mkdir -p $ltmpout/log/gsifit
tar czf $ltmpout/log/gsifit/gsifit_${Eyyyymmddhh}.tar.gz fort.*
EOF
fi
cat >> cycle_52.sh << EOF
rm -f $ltmpobs/obs_input.*
mv -f obs_input.* $ltmpobs
EOF

$MPIBIN/mpiexec -host ${node_m[$mmean]} bash cycle_52.sh &
wait

#-------------------------------------------------------------------------------
fi
#-------------------------------------------------------------------------------

if [ "$SHAREDISK" = '0' ] && [ "$nnodes" -gt 1 ]; then
  echo
  echo "  Distribute files from node ${node_m[$mmean]} to other nodes..."

  cat $tmpnode/machinefile.node | grep -xv ${node_m[$mmean]} > cycle_53.machinefile
  if [ "$THIN_OPT" = '2' ]; then
    obsfilename='obs*.dat'
  elif [ "$OBSOPE_OPT" = '2' ]; then
    obsfilename='obs_input.*'
  fi

#------ Method 1: use $RCP directly copy files from node to node
#  cat > cycle_53.sh << EOF
#cd $ltmpobs
#rm -f ${obsfilename}
#$SCP ${node_m[$mmean]}:${ltmpobs}/${obsfilename} .
#EOF
#------ Method 2: copy files to host node first, then distribute
  mkdir -p $tmpdata
  rm -fr $tmpdata/*
  $MPIBIN/mpiexec -host ${node_m[$mmean]} \
    bash -c "$CPCOMM ${ltmpobs}/${obsfilename} ${HOSTPREFIX}${tmpdata}" &
  wait
  cat > cycle_53.sh << EOF
cd $ltmpobs
rm -f ${obsfilename}
$CPCOMM ${HOSTPREFIX}${tmpdata}/$obsfilename .
EOF
#------

  np=`cat cycle_53.machinefile | wc -l`
  $MPIBIN/mpiexec -machinefile cycle_53.machinefile -n $np bash cycle_53.sh &
  wait
  echo "  Done!"
fi

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
fi
#===============================================================================
# 6. Run obs operator for each member (using GSI)
p=6
if [ "${run_part[$p]}" = '1' ]; then
if [ "$OBSOPE_OPT" = '2' ]; then
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
for m in `seq $MEMBER`; do
  $MPIBIN/mpiexec -host ${node_m[$m]} bash run_gsi.sh \
                  $ltmpgsi/${name_m[$m]} $ltmpprog $ltmpfixgsi $ltmpfixcrtm \
                  $ltmpgfs/${name_m[$m]} $ltmpobs 2 $gsiwindow $THIN_OPT &
  sleep 0.05s
done
wait

#-------------------------------------------------------------------------------

if [ "$LOG_OPT" -le 2 ]; then
  mkdir -p $OUTDIR/log/gsi/${Syyyymmddhh}
fi
echo
ppnl=$((ppn*mem_nodes_gsi))
pcount=0
for m in `seq $MEMBER`; do
  np=`cat $tmpnode/machinefile.gsi.${name_m[$m]} | wc -l`
  if [ "${node_m[$m]}" = "${node[1]}" ]; then
    pcount=$((pcount+np))
    if [ "$pcount" -gt "$ppnl" ]; then
      echo "    wait..."
      wait
      pcount=$np
    fi
  fi
  echo "  run GSI for member ${name_m[$m]} on node ${nodes_gsi_m[$m]}"
  if [ "$LOG_OPT" -le 2 ]; then
    $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.gsi.${name_m[$m]} -n $np \
                    -wdir $ltmpgsi/${name_m[$m]} \
                    ./global_gsi > $OUTDIR/log/gsi/${Syyyymmddhh}/gsi_${name_m[$m]}.log 2>&1 &
  else
    $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.gsi.${name_m[$m]} -n $np \
                    -wdir $ltmpgsi/${name_m[$m]} \
                    ./global_gsi > /dev/null 2>&1 &
  fi
  sleep 0.05s
done
echo "    wait..."
wait

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
fi
#=================================================================================
# 7. (When use obsope) Run obs operator for each member
#    (When use GSI   ) Convert GSI diag files into LETKF obs2 format
p=7
if [ "${run_part[$p]}" = '1' ]; then
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI

#-------------------------------------------------------------------------------
if [ "$OBSOPE_OPT" = '1' ]; then
#-------------------------------------------------------------------------------

  cat > cycle_71.sh << EOF
mem="\$1"
cd $ltmpobsope/\$mem
rm -f obsin*.dat
for is in \$(seq $nslots); do
  isf="\$(printf '%02d' \$is)"
  ln -fs $ltmpobs/obs\${isf}.dat obsin\${isf}.dat
done
ln -fs $ltmpobs/obsmn.dat obsinmn.dat
ln -fs $ltmpprog/obsope .
if [ "\$mem" = "${name_m[$MEMBER]}" ]; then
  echo
  ./obsope
else
  ./obsope > /dev/null 2>&1
fi
EOF
  if [ "$OBSOUT_OPT" -le 2 ]; then
    cat >> cycle_71.sh << EOF
mkdir -p $ltmpout/obs/obs2/${Eyyyymmddhh}
cp -f obsout.dat $ltmpout/obs/obs2/${Eyyyymmddhh}/obs2_\${mem}.dat
EOF
  fi
  cat >> cycle_71.sh << EOF
mv -f obsout.dat $ltmpletkf/obs\${mem}.dat
EOF

#-------------------------------------------------------------------------------
elif [ "$OBSOPE_OPT" = '2' ]; then
#-------------------------------------------------------------------------------

  cat > cycle_71.sh << EOF
mem="\$1"
cd $ltmpgsi/\$mem
mkdir -p $ltmpout/obs/gsidiag/${Eyyyymmddhh}
listall="$gsi_obstypes"
for type in \$listall; do
  count=\$(ls pe*.\${type}_01* 2> /dev/null | wc -l)
  if [ "\$count" -gt 0 ]; then
    cat pe*.\${type}_01* > diag_\${type}
EOF
  if [ "$OBSOUT_OPT" -le 1 ]; then
    cat >> cycle_71.sh << EOF
    gzip -c diag_\${type} > $ltmpout/obs/gsidiag/${Eyyyymmddhh}/diag_\${type}_\${mem}.gz
EOF
  fi
  cat >> cycle_71.sh << EOF
  fi
done

mkdir -p $ltmpreaddiag/\$mem
rm -f $ltmpreaddiag/\${mem}/*
cd $ltmpreaddiag/\$mem
ln -s $ltmpprog/readdiag_conv .
ln -fs $ltmpgsi/\${mem}/diag_\${type} obsin.dat
EOF
  if [ "$LOG_OPT" -le 2 ]; then
    cat >> cycle_71.sh << EOF
mkdir -p $ltmpout/log/readdiag/${Eyyyymmddhh}
./readdiag_conv > $ltmpout/log/readdiag/${Eyyyymmddhh}/readdiag_conv_\${mem}.log 2>&1
EOF
  else
    cat >> cycle_71.sh << EOF
./readdiag_conv > /dev/null 2>&1
EOF
  fi
  if [ "$OBSOUT_OPT" -le 2 ]; then
    cat >> cycle_71.sh << EOF
mkdir -p $ltmpout/obs/obs2/${Eyyyymmddhh}
cp -f obsout.dat $ltmpout/obs/obs2/${Eyyyymmddhh}/obs2_\${mem}.dat
EOF
  fi
  cat >> cycle_71.sh << EOF
mv -f obsout.dat $ltmpletkf/obs\${mem}.dat
EOF

#-------------------------------------------------------------------------------
fi
#-------------------------------------------------------------------------------

echo
ppnl=$ppn
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
  $MPIBIN/mpiexec -host ${node_m[$m]} bash cycle_71.sh ${name_m[$m]} &
  sleep 0.05s
done
echo "    wait..."
wait

#-------------------------------------------------------------------------------
# Modified LETKF program, do not need this distribution!!

#if [ "$SHAREDISK" = '0' ] && [ "$nnodes" -gt 1 ]; then
#  echo
#  echo "  Distribute files from node ${node_m[1]} to other nodes..."

#  cat $tmpnode/machinefile.node | grep -xv ${node_m[1]} > cycle_72.machinefile

##------ Method 1: use $RCP directly copy files from node to node
##  cat > cycle_72.sh << EOF
##cd $ltmpletkf
##if [ ! -s obs${name_m[1]}.dat ]; then
##  obsfilefind="\$(ls obs???.dat 2> /dev/null | head -n 1)"
##  if [ -n "\$obsfilefind" ]; then
##    ln -fs \$obsfilefind obs${name_m[1]}.dat
##  else
##    $SCP ${node_m[1]}:${$ltmpletkf}/obs${name_m[1]}.dat .
##  fi
##fi
##EOF
##------ Method 2: copy files to host node first, then distribute
#  mkdir -p $tmpdata
#  rm -fr $tmpdata/*
#  if [ "$nnodes" -gt "$MEMBER" ]; then # This criterion only works for cyclic member distribution
#    $MPIBIN/mpiexec -host ${node_m[1]} \
#      bash -c "$CPCOMM ${ltmpletkf}/obs${name_m[1]}.dat ${HOSTPREFIX}${tmpdata}" &
#    wait
#  fi
#  cat > cycle_72.sh << EOF
#cd $ltmpletkf
#if [ ! -s obs${name_m[1]}.dat ]; then
#  obsfilefind="\$(ls obs???.dat 2> /dev/null | head -n 1)"
#  if [ -n "\$obsfilefind" ]; then
#    ln -fs \$obsfilefind obs${name_m[1]}.dat
#  else
#    $CPCOMM ${HOSTPREFIX}${tmpdata}/obs${name_m[1]}.dat .
#  fi
#fi
#EOF
##------

#  np=`cat cycle_72.machinefile | wc -l`
#  $MPIBIN/mpiexec -machinefile cycle_72.machinefile -n $np bash cycle_72.sh &
#  wait
#  echo "  Done!"
#fi

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
#===============================================================================
# 8. LETKF
p=8
if [ "${run_part[$p]}" = '1' ]; then
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI

if [ "$SHAREDISK" = '0' ]; then
  cat > cycle_81.sh << EOF
cd $ltmpletkf
rm -f gues_me.grd gues_sp.grd infl_mul.grd
rm -f NOUT-*
ln -fs $ltmpprog/$letkfexec .
if [ -s "$ltmpobs/cdfm.grd" ]; then
  ln -fs $ltmpobs/cdfm.grd .
fi
if [ -s "$ltmpobs/cdfo.grd" ]; then
  ln -fs $ltmpobs/cdfo.grd .
fi
if [ -s "$ltmpobs/ppmask.grd" ]; then
  ln -fs $ltmpobs/ppmask.grd .
fi
EOF
  if [ "$ADAPTINFL" -ge 1 ]; then
    cat >> cycle_81.sh << EOF
prefixlen="\$(expr length \$(hostname))"
nodeinfl="${node_m[$mmean]}"
if [ "\$(hostname)" = "\${nodeinfl:0:\$prefixlen}" ]; then
  if [ -s "$ltmpout/infl/$INFLyyyymmddhh.grd" ]; then
    cp -f $ltmpout/infl/$INFLyyyymmddhh.grd infl_mul.grd
  fi
fi
EOF
  fi
  $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes bash cycle_81.sh &
  wait
else
  cd $ltmpletkf
  rm -f gues_me.grd gues_sp.grd infl_mul.grd
  rm -f NOUT-*
  ln -fs $ltmpprog/$letkfexec .
  if [ -s "$ltmpobs/cdfm.grd" ]; then
    ln -fs $ltmpobs/cdfm.grd .
  fi
  if [ -s "$ltmpobs/cdfo.grd" ]; then
    ln -fs $ltmpobs/cdfo.grd .
  fi
  if [ -s "$ltmpobs/ppmask.grd" ]; then
    ln -fs $ltmpobs/ppmask.grd .
  fi
  if [ "$ADAPTINFL" -ge 1 ] && [ -s "$ltmpout/infl/$INFLyyyymmddhh.grd" ]; then
    cp -f $ltmpout/infl/$INFLyyyymmddhh.grd infl_mul.grd
  fi
  cd $TMPMPI
fi

#-------------------------------------------------------------------------------

np=`cat $tmpnode/machinefile | wc -l`
echo
echo "  Run LETKF on $np processor cores..."

$MPIBIN/mpiexec -machinefile $tmpnode/machinefile -n $np \
                -wdir $ltmpletkf ./$letkfexec > /dev/null 2>&1 &
wait
echo "  Done!"

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
#===============================================================================
# 9. Convert LETKF output into GFS sig/sfc format
p=9
if [ "${run_part[$p]}" = '1' ]; then
part_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI

cat > cycle_91.sh << EOF
mem="\$1"
mkdir -p $ltmpssio/\$mem
cd $ltmpssio/\$mem
rm -f fort.*
ln -fs $ltmpprog/ss2grd .
ln -fs $ltmpprog/ss2grdp .
ln -fs $ltmpprog/grd2ss .
ln -fs $ltmpgfs/\${mem}/SIG.F${fhh_sl[$baseslot]} fort.11
ln -fs $ltmpgfs/\${mem}/SFC.F${fhh_sl[$baseslot]} fort.12
cp -f $ltmpanal/Eanal.sig fort.21
cp -f $ltmpanal/Eanal.sfc fort.22
ln -fs $ltmpletkf/anal\${mem}.grd fort.41
./grd2ss
mkdir -p $ltmpout/anal/\${mem}
mv -f fort.21 $ltmpout/anal/\${mem}/${Eyyyymmddhh}.sig
mv -f fort.22 $ltmpout/anal/\${mem}/${Eyyyymmddhh}.sfc
EOF
if [ "$OUT_OPT" -le 2 ]; then
  cat >> cycle_91.sh << EOF
mkdir -p $ltmpout/analg/\${mem}
mv -f $ltmpletkf/anal\${mem}.grd $ltmpout/analg/\${mem}/${Eyyyymmddhh}.grd
EOF
fi
if [ "$OUT_OPT" -le 1 ]; then
  cat >> cycle_91.sh << EOF
rm -f fort.*
ln -fs $ltmpout/anal/\${mem}/${Eyyyymmddhh}.sig fort.11
ln -fs $ltmpout/anal/\${mem}/${Eyyyymmddhh}.sfc fort.12
./ss2grdp
mkdir -p $ltmpout/analgp/\${mem}
mv -f fort.31 $ltmpout/analgp/\${mem}/${Eyyyymmddhh}.grd
ln -fs $ltmpout/gues/\${mem}/${Eyyyymmddhh}.sig fort.11
ln -fs $ltmpout/gues/\${mem}/${Eyyyymmddhh}.sfc fort.12
./ss2grdp
mkdir -p $ltmpout/guesgp/\${mem}
mv -f fort.31 $ltmpout/guesgp/\${mem}/${Eyyyymmddhh}.grd
EOF
fi
if [ "$OBSOUT_DIAG" = '1' ]; then
  cat >> cycle_91.sh << EOF
mkdir -p $ltmpout/obs/obsdiag/${Eyyyymmddhh}
if [ -e "$ltmpletkf/obsgues\${mem}.dat" ]; then
  mv -f $ltmpletkf/obsgues\${mem}.dat $ltmpout/obs/obsdiag/${Eyyyymmddhh}/obsgues_\${mem}.dat
fi
if [ -e "$ltmpletkf/obsanal\${mem}.dat" ]; then
  mv -f $ltmpletkf/obsanal\${mem}.dat $ltmpout/obs/obsdiag/${Eyyyymmddhh}/obsanal_\${mem}.dat
fi
EOF
fi

#-------------------------------------------------------------------------------

cat > cycle_91_node1.sh << EOF
mkdir -p $ltmpssio/mean
cd $ltmpssio/mean
rm -f fort.*
ln -fs $ltmpprog/ss2grd .
ln -fs $ltmpprog/ss2grdp .
ln -fs $ltmpprog/grd2ss .
EOF
if [ "$OUT_OPT" -le 6 ]; then
  cat >> cycle_91_node1.sh << EOF
ln -fs $ltmpgfs/mean/SIG.F${fhh_sl[$baseslot]} fort.11
ln -fs $ltmpgfs/mean/SFC.F${fhh_sl[$baseslot]} fort.12
cp -f $ltmpanal/Eanal.sig fort.21
cp -f $ltmpanal/Eanal.sfc fort.22
ln -fs $ltmpletkf/anal_me.grd fort.41
./grd2ss
mkdir -p $ltmpout/anal/mean
mv -f fort.21 $ltmpout/anal/mean/${Eyyyymmddhh}.sig
mv -f fort.22 $ltmpout/anal/mean/${Eyyyymmddhh}.sfc
EOF
fi
if [ "$OUT_OPT" -le 6 ]; then
  cat >> cycle_91_node1.sh << EOF
mkdir -p $ltmpout/analg/mean
mv -f $ltmpletkf/anal_me.grd $ltmpout/analg/mean/${Eyyyymmddhh}.grd
mkdir -p $ltmpout/analg/sprd
mv -f $ltmpletkf/anal_sp.grd $ltmpout/analg/sprd/${Eyyyymmddhh}.grd
EOF
fi
if [ "$OUT_OPT" -le 5 ]; then
  cat >> cycle_91_node1.sh << EOF
mkdir -p $ltmpout/guesg/mean
mv -f $ltmpletkf/gues_me.grd $ltmpout/guesg/mean/${Eyyyymmddhh}.grd
mkdir -p $ltmpout/guesg/sprd
mv -f $ltmpletkf/gues_sp.grd $ltmpout/guesg/sprd/${Eyyyymmddhh}.grd
EOF
fi
if [ "$OUT_OPT" -le 6 ]; then
  cat >> cycle_91_node1.sh << EOF
rm -f fort.*
ln -fs $ltmpout/anal/mean/${Eyyyymmddhh}.sig fort.11
ln -fs $ltmpout/anal/mean/${Eyyyymmddhh}.sfc fort.12
./ss2grdp
mkdir -p $ltmpout/analgp/mean
mv -f fort.31 $ltmpout/analgp/mean/${Eyyyymmddhh}.grd
EOF
fi
if [ "$OUT_OPT" -le 5 ]; then
  cat >> cycle_91_node1.sh << EOF
ln -fs $ltmpout/gues/mean/${Eyyyymmddhh}.sig fort.11
ln -fs $ltmpout/gues/mean/${Eyyyymmddhh}.sfc fort.12
./ss2grdp
mkdir -p $ltmpout/guesgp/mean
mv -f fort.31 $ltmpout/guesgp/mean/${Eyyyymmddhh}.grd
EOF
fi
if [ "$LOG_OPT" -le 3 ]; then
  cat >> cycle_91_node1.sh << EOF
mkdir -p $ltmpout/log/letkf
mv -f $ltmpletkf/NOUT-000 $ltmpout/log/letkf/letkf_${Eyyyymmddhh}.log
EOF
fi
if [ "$OBSOUT_OPT" -le 4 ]; then
  cat >> cycle_91_node1.sh << EOF
mkdir -p $ltmpout/obs/obsdep
mv -f $ltmpletkf/omb.dat $ltmpout/obs/obsdep/omb_$Eyyyymmddhh.dat
mv -f $ltmpletkf/oma.dat $ltmpout/obs/obsdep/oma_$Eyyyymmddhh.dat
EOF
fi
if [ "$OBSOUT_DIAG" = '1' ]; then
  cat >> cycle_91_node1.sh << EOF
mkdir -p $ltmpout/obs/obsdiag/${Eyyyymmddhh}
if [ -e "$ltmpletkf/obsgues_me.dat" ]; then
  mv -f $ltmpletkf/obsgues_me.dat $ltmpout/obs/obsdiag/${Eyyyymmddhh}/obsgues_mean.dat
fi
if [ -e "$ltmpletkf/obsanal_me.dat" ]; then
  mv -f $ltmpletkf/obsanal_me.dat $ltmpout/obs/obsdiag/${Eyyyymmddhh}/obsanal_mean.dat
fi
EOF
fi
if [ "$ADAPTINFL" -ge 1 ]; then
  cat >> cycle_91_node1.sh << EOF
mkdir -p $ltmpout/infl
cp -f $ltmpletkf/infl_mul.grd $ltmpout/infl/$EEyyyymmddhh.grd
EOF
fi

#-------------------------------------------------------------------------------

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
  $MPIBIN/mpiexec -host ${node_m[$m]} bash cycle_91.sh ${name_m[$m]} &
  sleep 0.05s
done

pcount=$((pcount+np))
if [ "$pcount" -gt "$ppnl" ]; then
  echo "    wait..."
  wait
fi
echo "  mean/sprd  on node ${node_m[$mmean]}"
$MPIBIN/mpiexec -host ${node_m[$mmean]} bash cycle_91_node1.sh &
echo "    wait..."
wait

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
#===============================================================================
# 10. Collect outputs (stageout)
p=10
if [ "${run_part[$p]}" = '1' ]; then
part_timing $p
#if [ "$MCYCLE" = '0' ] || [ "$MCYCLE" = '3' ]; then
#  echo "[$(now_format)] ${partname[$p]}" 1>&2
#  echo
#  printf " %2d. %-55s\n" ${p} "${partname[$p]}"
#else
#  echo "[$(now_format)] ${partname[$p]} <send to background>" 1>&2
#  echo
#  printf " %2d. %-55s\n" ${p} "${partname[$p]} <send to background>"
#fi
#-------------------------------------------------------------------------------
if [ "$SHAREDISK" = '0' ]; then
#-------------------------------------------------------------------------------

cd $TMPMPI
mkdir -p $tmpstageout
rm -f $tmpstageout/*

if [ "$OUT_OPT" -le 4 ]; then
  for m in `seq $MEMBER`; do
    echo "rm|anal/${name_m[$m]}/${Syyyymmddhh}.sig" >> $tmpstageout/out.${node_m[$m]}
    echo "rm|anal/${name_m[$m]}/${Syyyymmddhh}.sfc" >> $tmpstageout/out.${node_m[$m]}
    echo "cp|anal/${name_m[$m]}/${Eyyyymmddhh}.sig" >> $tmpstageout/out.${node_m[$m]}
    echo "cp|anal/${name_m[$m]}/${Eyyyymmddhh}.sfc" >> $tmpstageout/out.${node_m[$m]}
  done
fi
echo "mv|anal/mean/${Eyyyymmddhh}.sig" >> $tmpstageout/out.${node_m[$mmean]}
echo "mv|anal/mean/${Eyyyymmddhh}.sfc" >> $tmpstageout/out.${node_m[$mmean]}
for m in `seq $msprd`; do
  echo "mv|analg/${name_m[$m]}/${Eyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$m]}
done
for m in `seq $mmean`; do
  echo "mv|analgp/${name_m[$m]}/${Eyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$m]}
done
for m in `seq $mmean`; do
  echo "mv|gues/${name_m[$m]}/${Eyyyymmddhh}.sig" >> $tmpstageout/out.${node_m[$m]}
  echo "mv|gues/${name_m[$m]}/${Eyyyymmddhh}.sfc" >> $tmpstageout/out.${node_m[$m]}
done
for m in `seq $msprd`; do
  echo "mv|guesg/${name_m[$m]}/${Eyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$m]}
done
for m in `seq $mmean`; do
  echo "mv|guesgp/${name_m[$m]}/${Eyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$m]}
done
if [ "$ADAPTINFL" -ge 1 ]; then
  echo "rm|infl/${INFLyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mmean]}
  echo "cp|infl/${EEyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mmean]}
fi

cat >> $tmpstageout/out.${node_m[$mmean]} << EOF
mv|obs/obsinput/${Eyyyymmddhh}/obs_input.tar.gz
mv|obs/gsidiag/${Eyyyymmddhh}/diag_conv_${name_m[$mmean]}.gz
mv|obs/obsdep/omb_${Eyyyymmddhh}.dat
mv|obs/obsdep/oma_${Eyyyymmddhh}.dat
mv|log/superob/superob_${Eyyyymmddhh}.log
mv|log/superob_d/superob_${Eyyyymmddhh}.tar.gz
mv|log/gsifit/gsifit_${Eyyyymmddhh}.tar.gz
mv|log/letkf/letkf_${Eyyyymmddhh}.log
EOF
if [ "$THIN_OPT" = '2' ]; then
  for is in `seq $nslots`; do
    echo "mv|obs/superob/obs${Eyyyymmddhh}/${obsfile_sl[$is]}" >> $tmpstageout/out.${node_m[$mmean]}
  done
fi
if [ "$LOG_OPT" -le 2 ]; then
  for m in `seq $mmean`; do
    echo "mv|log/readdiag/${Eyyyymmddhh}/readdiag_conv_${name_m[$m]}.log" >> $tmpstageout/out.${node_m[$m]}
  done
fi
if [ "$OBSOUT_OPT" -le 1 ]; then
  for m in `seq $mmean`; do
    echo "mv|obs/gsidiag/${Eyyyymmddhh}/diag_conv_${name_m[$m]}.gz" >> $tmpstageout/out.${node_m[$m]}
  done
fi
if [ "$OBSOUT_OPT" -le 2 ]; then
  for m in `seq $mmean`; do
    echo "mv|obs/obs2/${Eyyyymmddhh}/obs2_${name_m[$m]}.dat" >> $tmpstageout/out.${node_m[$m]}
  done
fi
if [ "$OBSOUT_DIAG" = '1' ]; then
  for m in `seq $mmean`; do
    echo "mv|obs/obsdiag/${Eyyyymmddhh}/obsgues_${name_m[$m]}.dat" >> $tmpstageout/out.${node_m[$m]}
    echo "mv|obs/obsdiag/${Eyyyymmddhh}/obsanal_${name_m[$m]}.dat" >> $tmpstageout/out.${node_m[$m]}
  done
fi

#-------------------------------------------------------------------------------

#if [ -e "stageout_bg.running" ]; then
#  echo
#  echo "  wait for previous stageout to finish..."
#fi
#while [ -e "stageout_bg.running" ]; do
#  sleep 1s
#done

if [ "$MCYCLE" = '0' ] || [ "$MCYCLE" = '3' ]; then # clean stageout
  stageout $ltmpout 0
  $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes \
                  rm -fr $LTMP1/${tmpsubdir} $LTMP2/${tmpsubdir} &
  wait
else # not clean stageout, keep local files for the next cycle
  stageout $ltmpout 0 &
  wait
#  stageout $ltmpout 1
fi

#-------------------------------------------------------------------------------
elif [ "$SHAREDISK" = '1' ]; then
#-------------------------------------------------------------------------------

if [ "$MCYCLE" = '0' ] || [ "$MCYCLE" = '3' ]; then # clean stageout
  rm -fr $ltmprun1 $ltmprun2
fi

#-------------------------------------------------------------------------------
fi
#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
#===============================================================================

#rm -fr $TMPMPI

echo
echo " +----------------------------------------------------------------+"
echo " |                GFS-LETKF successfully completed                |"
echo " +----------------------------------------------------------------+"
echo

echo "[$(now_format)] Finish cycle.sh $@" 1>&2
exit 0

#===============================================================================
