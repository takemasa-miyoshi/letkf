#!/bin/bash
#===============================================================================
#
#  Script to post-process the obsope outputs.
#  December 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 7)); then
  cat >&2 << EOF

[post_obsope.sh]

Usage: $0 MYRANK STIME ATIME MEM TMPDIR LOG_OPT OUT_OPT

  MYRANK  My rank number
  STIME
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  MEM     Name of the ensemble member
  TMPDIR  Temporary directory to run the program
  LOG_OPT
  OUT_OPT

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
ATIME="$1"; shift
MEM="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"; shift
OUT_OPT="$1"

#===============================================================================

if ((LOG_OPT <= 4 && MYRANK == 0)); then
  if [ -f "$TMPDIR/obsope.conf" ]; then
    mv -f $TMPDIR/obsope.conf $TMPOUT/${ATIME}/log/obsope/obsope.conf
    mv -f $TMPDIR/LOG.pe000000 $TMPOUT/${ATIME}/log/obsope/LOG.pe000000
  fi
fi

if ((OUT_OPT >= 3)); then
  if [ -d "$TMPOUT/${STIME}/hist/${MEM}" ]; then
    rm -f $TMPOUT/${STIME}/hist/${MEM}/*
  fi
fi

#===============================================================================

exit 0
