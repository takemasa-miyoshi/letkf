#!/bin/bash
#===============================================================================
#
#  Date and time utility (using unix 'date' command)
#  May 2014, Guo-Yuan Lien
#
#===============================================================================

function datetime {
#-------------------------------------------------------------------------------
# Return a full timestamp from a reference (partial) timestamp
# and an optional time increment.
#
# Usage: datetime YYYY[MMDDHHIISS] [INC] [UNIT]
#
#   YYYYMMDDHHIISS  Reference (partial) timestamp
#                     YYYY: year
#                     MM:   month
#                     DD:   day
#                     HH:   hour
#                     II:   minute
#                     SS:   second
#   INC             Increment
#                   (default: 0)
#   UNIT            Unit of the increment
#                     y: year
#                     m: month
#                     d: day
#                     h: hour
#                     i: minute
#                     s: second
#                   (default: s)
#
# Return:
#   YYYYMMDDHHIISS  Resultant full timestamp
#-------------------------------------------------------------------------------

if [ "$#" -lt 1 ]; then
  echo "[Error] function 'datetime': Insufficient arguments." 1>&2
  exit 1
fi

local YYYYMMDDHHIISS="$1"
local INC="${2:-0}"
local UNIT="${3:-second}"

datecmd="/bin/date"

#-------------------------------------------------------------------------------

if [ "$UNIT" = 'y' ]; then
  UNIT='year'
fi
if [ "$UNIT" = 'm' ]; then
  UNIT='month'
fi
if [ "$UNIT" = 'd' ]; then
  UNIT='day'
fi
if [ "$UNIT" = 'h' ]; then
  UNIT='hour'
fi
if [ "$UNIT" = 'i' ]; then
  UNIT='minute'
fi
if [ "$UNIT" = 's' ]; then
  UNIT='second'
fi

local YYYY="${YYYYMMDDHHIISS:0:4}"
local MM="${YYYYMMDDHHIISS:4:2}"
local DD="${YYYYMMDDHHIISS:6:2}"
local HH="${YYYYMMDDHHIISS:8:2}"
local II="${YYYYMMDDHHIISS:10:2}"
local SS="${YYYYMMDDHHIISS:12:2}"

MM="${MM:-01}"
DD="${DD:-01}"
HH="${HH:-00}"
II="${II:-00}"
SS="${SS:-00}"

res=$($datecmd -ud "$INC $UNIT ${YYYY}-${MM}-${DD} ${HH}:${II}:${SS}" +'%Y%m%d%H%M%S')
if [ "$?" -eq 0 ]; then
  echo $res
else
  echo "[Error] datetime: Wrong date calculation." 1>&2
  exit 1
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function datetimegrads {
#-------------------------------------------------------------------------------
# Format date and time for GrADS
#
# Usage: datetimegrads YYYY[MMDDHHIISS]
#
#   YYYYMMDDHHIISS  Reference (partial) timestamp
#                     YYYY: year
#                     MM:   month
#                     DD:   day
#                     HH:   hour
#                     II:   minute
#                     SS:   second
#-------------------------------------------------------------------------------

if [ "$#" -lt 1 ]; then
  echo "[Error] function 'datetimegrads': Insufficient arguments." 1>&2
  exit 1
fi

local YYYYMMDDHHIISS=$(datetime $1)

#-------------------------------------------------------------------------------

local YYYY="${YYYYMMDDHHIISS:0:4}"
local MM="${YYYYMMDDHHIISS:4:2}"
local DD="${YYYYMMDDHHIISS:6:2}"
local HH="${YYYYMMDDHHIISS:8:2}"
local II="${YYYYMMDDHHIISS:10:2}"
local mmm

if [ "$MM" = '01' ]; then
  mmm='jan'
elif [ "$MM" = '02' ]; then
  mmm='feb'
elif [ "$MM" = '03' ]; then
  mmm='mar'
elif [ "$MM" = '04' ]; then
  mmm='apr'
elif [ "$MM" = '05' ]; then
  mmm='may'
elif [ "$MM" = '06' ]; then
  mmm='jun'
elif [ "$MM" = '07' ]; then
  mmm='jul'
elif [ "$MM" = '08' ]; then
  mmm='aug'
elif [ "$MM" = '09' ]; then
  mmm='sep'
elif [ "$MM" = '10' ]; then
  mmm='oct'
elif [ "$MM" = '11' ]; then
  mmm='nov'
elif [ "$MM" = '12' ]; then
  mmm='dec'
fi

if [ "$II" = '00' ]; then
  echo "${HH}Z${DD}${mmm}${YYYY}"
else
  echo "${HH}:${II}Z${DD}${mmm}${YYYY}"
fi

#-------------------------------------------------------------------------------
}

#===============================================================================
