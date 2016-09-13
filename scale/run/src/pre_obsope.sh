#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of obsope run; for each member.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 3)); then
  cat >&2 << EOF

[pre_obsope.sh]

Usage: $0 MYRANK ATIME MEM

  MYRANK  My rank number (not used)
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  MEM     Name of the ensemble member

EOF
  exit 1
fi

MYRANK="$1"; shift
ATIME="$1"; shift
MEM="$1"

#===============================================================================

if [ "$MEM" != 'mean' ]; then
  mkdir -p $TMPOUT/${ATIME}/obsgues/${MEM}
fi

#===============================================================================

exit 0
