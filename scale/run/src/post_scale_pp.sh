#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 5)); then
  cat >&2 << EOF

[post_scale_pp.sh] Post-process the scale-rm_init outputs.

Usage: $0 MYRANK STIME MEM TMPDIR LOG_OPT [SCPCALL]

  MYRANK   My rank number (not used)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MEM
  TMPDIR   Temporary directory to run the model
  LOG_OPT
  SCPCALL

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
MEM="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"; shift
SCPCALL="${1:-cycle}"

#===============================================================================

if [ "$SCPCALL" = 'cycle' ]; then
  if ((LOG_OPT <= 2)); then
    if [ -f "$TMPDIR/pp.conf" ]; then
      mv -f $TMPDIR/pp.conf $TMPOUT/${STIME}/log/scale_pp/${MEM}_pp.conf
    fi
  fi
elif [ "$SCPCALL" = 'fcst' ]; then
  if ((LOG_OPT <= 2)); then
    if [ -f "$TMPDIR/pp.conf" ]; then
      mv -f $TMPDIR/pp.conf $TMPOUT/${STIME}/log/${SCPCALL}_scale_pp/${MEM}_pp.conf
    fi
  fi
fi

#===============================================================================

exit 0
