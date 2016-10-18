#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main
. src/func_datetime.sh

if (($# < 7)); then
  cat >&2 << EOF

[post_scale.sh] Post-process the SCALE model outputs.

Usage: $0 MYRANK STIME MEM FCSTLEN TMPDIR LOG_OPT OUT_OPT [SCPCALL] [DELETE_MEMBER]

  MYRANK   My rank number (not used)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MEM      Name of the ensemble member
  FCSTLEN  Forecast length (second)
  TMPDIR   Temporary directory to run the model
  LOG_OPT
  OUT_OPT
  SCPCALL  Called from which script? (fcst/cycle)
  DELETE_MEMBER

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
MEM="$1"; shift
FCSTLEN="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"; shift
OUT_OPT="$1"; shift
SCPCALL="${1:-cycle}"; shift
DELETE_MEMBER="${1:-0}"

ATIME=$(datetime $STIME $LCYCLE s)

restartbaselen=27  # 7 + 20

#===============================================================================

if [ "$SCPCALL" = 'cycle' ]; then

  MEMtmp=$MEM
  if [ "$MEM" = 'mean' ]; then
    MEMtmp='meanf'
  fi
  mkdir -p $TMPOUT/${STIME}/hist/${MEMtmp}
  mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/hist/${MEMtmp}
  mkdir -p $TMPOUT/${ATIME}/gues/${MEMtmp}
  file_prefix=$(cd $TMPDIR ; ls restart*.nc | head -n 1) # pick up the first restart output. ###### TO DO: explicitly calculate the time string???
  for ifile in $(cd $TMPDIR ; ls ${file_prefix:0:$restartbaselen}*.nc); do
    mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/gues/${MEMtmp}/init${ifile:$restartbaselen}
  done

  if ((LOG_OPT <= 3)); then
    if [ -f "$TMPDIR/run.conf" ]; then
      mv -f $TMPDIR/run.conf $TMPOUT/${STIME}/log/scale/${MEM}_run.conf
    fi
  fi

  if ((DELETE_MEMBER == 1)) && [ "$MEM" != 'mean' ]; then
    if [ -d "$TMPOUT/${STIME}/gues/${MEM}" ]; then
      rm -f $TMPOUT/${STIME}/gues/${MEM}/*
    fi
    if [ -d "$TMPOUT/${STIME}/anal/${MEM}" ]; then
      rm -f $TMPOUT/${STIME}/anal/${MEM}/*
    fi
  fi

  if ((MYRANK == 0)); then
    if [ -f "$TMPDIR/../latlon_domain_catalogue.txt" ]; then
      mv -f $TMPDIR/../latlon_domain_catalogue.txt $TMPOUT/${STIME}/log/scale/latlon_domain_catalogue.txt
    fi
  fi

elif [ "$SCPCALL" = 'fcst' ]; then

  mkdir -p $TMPOUT/${STIME}/fcst/${MEM}
  mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/fcst/${MEM}
  if ((OUT_OPT <= 1)); then
    file_prefix=$(cd $TMPDIR ; ls restart*.nc | head -n 1) # pick up the first restart output. ###### TO DO: explicitly calculate the time string???
    for ifile in $(cd $TMPDIR ; ls ${file_prefix:0:$restartbaselen}*.nc); do
      mv -f ${TMPDIR}/${ifile} $TMPOUT/${STIME}/fcst/${MEM}/init_$(datetime ${STIME} $FCSTLEN s)${ifile:$restartbaselen}
    done
  fi

  if ((LOG_OPT <= 3)); then
    if [ -f "$TMPDIR/run.conf" ]; then
      mv -f $TMPDIR/run.conf $TMPOUT/${STIME}/log/${SCPCALL}_scale/${MEM}_run.conf
    fi
  fi

  if ((MYRANK == 0)); then
    if [ -f "$TMPDIR/../latlon_domain_catalogue.txt" ]; then
      mv -f $TMPDIR/../latlon_domain_catalogue.txt $TMPOUT/${STIME}/log/${SCPCALL}_scale/latlon_domain_catalogue.txt
    fi
  fi

fi

#===============================================================================

exit 0
