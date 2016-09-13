#!/bin/bash
#===============================================================================
#
#  Date and time utilities
#  May 2014, Guo-Yuan Lien
#
#  *Require source 'config.main' first.
#
#===============================================================================

Unix_date_cmd=1  # 0: Use the built-in 'datetime' program
                 # 1: Use the Unix command 'date'

#===============================================================================

datetime () {
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
  echo "[Error] $FUNCNAME: Insufficient arguments." 1>&2
  exit 1
fi

local YYYYMMDDHHIISS="$1"
local INC="${2:-0}"
local UNIT="${3:-s}"

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

if [ "$Unix_date_cmd" -eq 0 ]; then

  if [ "$STG_TYPE" = 'K' ] || [ "$STG_TYPE" = 'K_rankdir' ]; then
    local datetime_prog="$TMPDAT_DIR/exec/datetime"
  else
    local datetime_prog="$COMMON_DIR/datetime"
  fi
  $datetime_prog << EOF
$YYYY
$MM
$DD
$HH
$II
$SS
$INC
$UNIT
EOF

elif [ "$Unix_date_cmd" -eq 1 ]; then

  if [ "$UNIT" = 'y' ]; then
    UNIT='year'
  elif [ "$UNIT" = 'm' ]; then
    UNIT='month'
  elif [ "$UNIT" = 'd' ]; then
    UNIT='day'
  elif [ "$UNIT" = 'h' ]; then
    UNIT='hour'
  elif [ "$UNIT" = 'i' ]; then
    UNIT='minute'
  elif [ "$UNIT" = 's' ]; then
    UNIT='second'
  fi
  date -ud "$INC $UNIT ${YYYY}-${MM}-${DD} ${HH}:${II}:${SS}" +'%Y%m%d%H%M%S'

fi

if [ $? -ne 0 ]; then
  local ierr=$?
  echo "[Error] $FUNCNAME: Wrong date calculation." 1>&2
  exit $ierr
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

datetime_fmt () {
#-------------------------------------------------------------------------------
# Print formatted date and time
#
# Usage: datetime_fmt YYYY[MMDDHHIISS]
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
  echo "[Error] $FUNCNAME: Insufficient arguments." 1>&2
  exit 1
fi

local YYYYMMDDHHIISS=$(datetime $1)

#-------------------------------------------------------------------------------

local YYYY="${YYYYMMDDHHIISS:0:4}"
local MM="${YYYYMMDDHHIISS:4:2}"
local DD="${YYYYMMDDHHIISS:6:2}"
local HH="${YYYYMMDDHHIISS:8:2}"
local II="${YYYYMMDDHHIISS:10:2}"
local SS="${YYYYMMDDHHIISS:12:2}"

echo "${YYYY}-${MM}-${DD} ${HH}:${II}:${SS}"

#-------------------------------------------------------------------------------
}

#===============================================================================

datetime_grads () {
#-------------------------------------------------------------------------------
# Print formatted date and time for GrADS
#
# Usage: datetime_grads YYYY[MMDDHHII]
#
#   YYYYMMDDHHII  Reference (partial) timestamp
#                   YYYY: year
#                   MM:   month
#                   DD:   day
#                   HH:   hour
#                   II:   minute
#-------------------------------------------------------------------------------

if [ "$#" -lt 1 ]; then
  echo "[Error] $FUNCNAME: Insufficient arguments." 1>&2
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

datetime_now () {
#-------------------------------------------------------------------------------
# Print current time
#
# Usage: datetime_now
#-------------------------------------------------------------------------------

  date +'%Y-%m-%d %H:%M:%S'

#-------------------------------------------------------------------------------
}

#===============================================================================
