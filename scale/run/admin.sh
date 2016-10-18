#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

. config.main
res=$? && ((res != 0)) && exit $res

#. src/func_datetime.sh

#-------------------------------------------------------------------------------

if (($# < 7)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

SCPNAME="$1"; shift
STIME="$1"; shift
ETIME="$1"; shift
TIME_DT="$1"; shift
TIME_DT_DYN="$1"; shift
NNODES="$1"; shift
WTIME_L="$1"

CONFIG='realtime_v160405_d1'

if [ "$ETIME" = '-' ]; then
  ETIME="$STIME"
fi

#-------------------------------------------------------------------------------

if [ "$SCPNAME" = 'cycle' ]; then
  DATA_BDY_WRF="ncepgfs_wrf_da"
else
  DATA_BDY_WRF="ncepgfs_wrf"
fi
#cat config/${CONFIG}/config.main.K_micro | \
cat config/${CONFIG}/config.main.K | \
    sed -e "s/<DATA_BDY_WRF>/${DATA_BDY_WRF}/g" | \
    sed -e "s/<NNODES>/${NNODES}/g" \
    > config.main

cat config/${CONFIG}/config.${SCPNAME} | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<ETIME>/${ETIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" \
    > config.${SCPNAME}

cat config/${CONFIG}/config.nml.scale | \
    sed -e "s/<TIME_DT>/${TIME_DT}/g" | \
    sed -e "s/<TIME_DT_DYN>/${TIME_DT_DYN}/g" \
    > config.nml.scale

#-------------------------------------------------------------------------------

#./${SCPNAME}_K_micro.sh > ${SCPNAME}_K.log 2>&1
./${SCPNAME}_K.sh > ${SCPNAME}_K.log 2>&1
res=$? && ((res != 0)) && exit $res

jobname="${SCPNAME}_${SYSNAME}"
jobid=$(grep 'pjsub Job' ${SCPNAME}_K.log | cut -d ' ' -f6)

#-------------------------------------------------------------------------------

if [ ! -s "${jobname}.o${jobid}" ] || [ ! -s "${jobname}.e${jobid}" ] || \
   [ ! -s "${jobname}.i${jobid}" ] || [ ! -s "${jobname}.s${jobid}" ]; then
  exit 101
elif [ -n "$(grep 'ERR.' ${jobname}.e${jobid})" ]; then
  exit 102
elif [ -n "$(grep 'terminated' ${jobname}.e${jobid})" ]; then
  exit 103
#elif [ ! -s "${jobname}.s${jobid}" ]; then
#  exit 104
#elif [ "$(tail -n 1 ${jobname}.s${jobid})" != "---(Stage-Out Error Information)---" ]; then
#  exit 105
fi

rm -f ${SCPNAME}_job.sh
rm -f ${jobname}.o${jobid}
rm -f ${jobname}.e${jobid}
rm -f ${jobname}.s${jobid}
rm -f ${jobname}.i${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_${SCPNAME}_${STIME} exp

#-------------------------------------------------------------------------------

exit 0
