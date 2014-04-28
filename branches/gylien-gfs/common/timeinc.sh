#!/bin/sh
set -e

function timeinc6hr
{
if test $# -ne 4
then
  echo "USAGE: $0 yyyy mm dd hh"
  return 1
fi

local YYYY
local MM
local DD
local HH
local ITMP

YYYY=$1
MM=$2
DD=$3
HH=$4

# Increment date
HH=`expr $HH + 6`
if test $HH -lt 10
then
  HH=0$HH
elif test $HH -gt 23
then
  HH=00
  DD=`expr $DD + 1`
  if test $DD -lt 10
  then
    DD=0$DD
  elif test $DD -eq 29
  then
    ITMP=`expr $YYYY % 4` || test 1 -eq 1
    if test $MM -eq 02 -a $ITMP -ne 0
    then
      DD=01
      MM=03
    fi
  elif test $DD -eq 30
  then
    ITMP=`expr $YYYY % 4` || test 1 -eq 1
    if test $MM -eq 02 -a $ITMP -eq 0
    then
      DD=01
      MM=03
    fi
  elif test $DD -eq 31
  then
    if test $MM -eq 04 -o $MM -eq 06
    then
      DD=01
      MM=`expr $MM + 1`
      MM=0$MM
    elif test $MM -eq 09 -o $MM -eq 11
    then
      DD=01
      MM=`expr $MM + 1`
    fi
  elif test $DD -eq 32
  then
    DD=01
    MM=`expr $MM + 1`
    if test $MM -lt 10
    then
      MM=0$MM
    fi
  fi
  if test $MM -gt 12
  then
    MM=01
    YYYY=`expr $YYYY + 1`
  fi
fi
#
# Outputs
#
echo $YYYY$MM$DD$HH
return 0
}

