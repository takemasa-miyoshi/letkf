#!/bin/bash
#===============================================================================

function ndarray_get {
#-------------------------------------------------------------------------------

local USAGE
read -d '' USAGE << EOF
'$FUNCNAME $@'
Usage: $FUNCNAME ARRAY INDEX1 INDEX2 ...
       The array dimensions are set in 'ARRAY_dimension' variable.
       ARRAY_dimension='DIM1 DIM2 ...'
  ARRAY   Name of the array
  INDEX1  Index of the first dimension
  INDEX2  Index of the second dimension
  ...
Return: Value of the array element
EOF

if [ -z "$1" ]; then
  echo "[Error] $FUNCNAME: The name of array is missing." 1>&2
#  exit 1
  return 1
fi

local array_dim_var="${1}_dimension"
local array_dim=(${!array_dim_var})
if [ -z "$array_dim" ]; then
  echo "[Error] $FUNCNAME: The dimensions of array '$ARRAY' have not been set." 1>&2
#  exit 1
  return 1
fi

local i
local ii
local idx
for ((i=0; i<${#array_dim[@]}; i++)); do
  ii=$((i + 2))
  if [ -z "${!ii}" ]; then
    echo "[Error] $FUNCNAME: DIM$((i+1)) is missing." 1>&2
    echo "$USAGE" 1>&2
#    exit 1
    return 1
  fi
  if ! [ "${!ii}" -eq "${!ii}" 2> /dev/null ]; then
    echo "[Error] $FUNCNAME: DIM$((i+1)) is not a number." 1>&2
    echo "$USAGE" 1>&2
#    exit 1
    return 1
  fi
  if ! [ "${array_dim[$i]}" -eq "${array_dim[$i]}" 2> /dev/null ]; then
    echo "[Error] $FUNCNAME: The dimensions of array '$ARRAY' are incorrect." 1>&2
    echo "$USAGE" 1>&2
#    exit 1
    return 1
  fi
  if [ "${!ii}" -lt 1 ] || [ "${!ii}" -gt "${array_dim[$i]}" ]; then
    echo "[Error] $FUNCNAME: DIM$((i+1)) should be within [1, ${array_dim[$i]}]." 1>&2
    echo "$USAGE" 1>&2
#    exit 1
    return 1
  fi

  if [ "$i" -eq 0 ]; then
    idx=$((${!ii} - 1))
  else
    idx=$((idx * ${array_dim[$i]} + ${!ii} - 1))
  fi
done

local expression="$1[$idx]"
echo "${!expression}"

#-------------------------------------------------------------------------------
}

#===============================================================================

function ndarray_set {
#-------------------------------------------------------------------------------

local USAGE
read -d '' USAGE << EOF
'$FUNCNAME $@'
Usage: $FUNCNAME ARRAY INDEX1 INDEX2 ... VALUE
       The array dimensions are set in 'ARRAY_dimension' variable.
       ARRAY_dimension='DIM1 DIM2 ...'
  ARRAY   Name of the array
  INDEX1  Index of the first dimension
  INDEX2  Index of the second dimension
  ...
  VALUE   Value of the array element
EOF

if [ -z "$1" ]; then
  echo "[Error] $FUNCNAME: The name of array is missing." 1>&2
#  exit 1
  return 1
fi

local array_dim_var="${1}_dimension"
local array_dim=(${!array_dim_var})
if [ -z "$array_dim" ]; then
  echo "[Error] $FUNCNAME: The dimensions of array '$ARRAY' have not been set." 1>&2
#  exit 1
  return 1
fi

local i
local ii
local idx
for ((i=0; i<${#array_dim[@]}; i++)); do
  ii=$((i + 2))
  if [ -z "${!ii}" ]; then
    echo "[Error] $FUNCNAME: DIM$((i+1)) is missing." 1>&2
    echo "$USAGE" 1>&2
#    exit 1
    return 1
  fi
  if ! [ "${!ii}" -eq "${!ii}" 2> /dev/null ]; then
    echo "[Error] $FUNCNAME: DIM$((i+1)) is not a number." 1>&2
    echo "$USAGE" 1>&2
#    exit 1
    return 1
  fi
  if ! [ "${array_dim[$i]}" -eq "${array_dim[$i]}" 2> /dev/null ]; then
    echo "[Error] $FUNCNAME: The dimensions of array '$ARRAY' are incorrect." 1>&2
    echo "$USAGE" 1>&2
#    exit 1
    return 1
  fi
  if [ "${!ii}" -lt 1 ] || [ "${!ii}" -gt "${array_dim[$i]}" ]; then
    echo "[Error] $FUNCNAME: DIM$((i+1)) should be within [1, ${array_dim[$i]}]." 1>&2
    echo "$USAGE" 1>&2
#    exit 1
    return 1
  fi

  if [ "$i" -eq 0 ]; then
    idx=$((${!ii} - 1))
  else
    idx=$((idx * ${array_dim[$i]} + ${!ii} - 1))
  fi
done

ii=$((${#array_dim[@]} + 2))
local expression="$1[$idx]"
eval $expression=${!ii}

#-------------------------------------------------------------------------------
}

#===============================================================================
