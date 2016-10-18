#!/bin/bash
#===============================================================================
#
#  Common utilities (using built-in 'datetime' program)
#  August 2014, Guo-Yuan Lien
#
#  *Require source 'config.main' first.
#
#===============================================================================

safe_init_tmpdir () {
#-------------------------------------------------------------------------------
# Safely initialize a temporary directory
#
# Usage: safe_init_tmpdir DIRNAME
#
#   DIRNAME  The temporary directory
#-------------------------------------------------------------------------------

local DIRNAME="$1"

#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Warning] $FUNCNAME: '\$DIRNAME' is not set." >&2
  exit 1
fi

mkdir -p $DIRNAME || exit $?

if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." >&2
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." >&2
  exit 1
fi

rm -fr $DIRNAME/* || exit $?

#-------------------------------------------------------------------------------
}

#===============================================================================

safe_rm_tmpdir () {
#-------------------------------------------------------------------------------
# Safely remove a temporary directory
#
# Usage: safe_rm_tmpdir DIRNAME
#
#   DIRNAME  The temporary directory
#-------------------------------------------------------------------------------

local DIRNAME="$1"

#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '\$DIRNAME' is not set." >&2
  exit 1
fi
if [ ! -e "$DIRNAME" ]; then
  return 0
fi
if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." >&2
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." >&2
  exit 1
fi

rm -fr $DIRNAME
res=$? && ((res != 0)) && exit $res

#-------------------------------------------------------------------------------
}

#===============================================================================

rev_path () {
#-------------------------------------------------------------------------------
# Compose the reverse path of a path
#
# Usage: rev_path PATH
#
#   PATH  The forward path
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local path="$1"

#-------------------------------------------------------------------------------

local rpath='.'
local base
while [ "$path" != '.' ]; do
  base=$(basename $path)
  res=$? && ((res != 0)) && exit $res
  path=$(dirname $path)
  if [ "$base" = '..' ]; then
    if [ -d "$path" ]; then
      rpath="$rpath/$(basename $(cd $path && pwd))"
    else
      echo "[Error] $FUNCNAME: Error in reverse path search." 1>&2
      exit 1
    fi
  elif [ "$base" != '.' ]; then
    rpath="$rpath/.."
  fi
done
if [ ${rpath:0:2} = './' ]; then
  echo ${rpath:2}
else
  echo $rpath
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

mpirunf () {
#-------------------------------------------------------------------------------
# Submit a MPI job according to nodefile
#
# Usage: mpirunf NODEFILE PROG [ARGS]
#
#   NODEFILE  Name of nodefile (omit the directory $NODEFILE_DIR)
#   PROG      Program
#   ARGS      Arguments passed into the program
#
# Other input variables:
#   $NODEFILE_DIR  Directory of nodefiles
#-------------------------------------------------------------------------------

if (($# < 2)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local NODEFILE="$1"; shift
local PROG="$1"; shift
local CONF="$1"; shift
local STDOUT="$1"; shift
local ARGS="$@"

progbase=$(basename $PROG)
progdir=$(dirname $PROG)

#-------------------------------------------------------------------------------

if [ "$MPI_TYPE" = 'sgimpt' ]; then

  local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

  $MPIRUN -d $progdir $HOSTLIST 1 ./$progbase $CONF $STDOUT $ARGS
#  $MPIRUN -d $progdir $HOSTLIST 1 omplace -nt ${THREADS} ./$progbase $CONF $STDOUT $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $MPIRUN -d $progdir $HOSTLIST 1 ./$progbase $CONF $STDOUT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'openmpi' ]; then

  NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)

  $MPIRUN -np $NNP -hostfile ${NODEFILE_DIR}/${NODEFILE} -wdir $progdir ./$progbase $CONF $STDOUT $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $$MPIRUN -np $NNP -hostfile ${NODEFILE_DIR}/${NODEFILE} -wdir $progdir ./$progbase $CONF $STDOUT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'K' ]; then

  NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)

  if [ "$STG_TYPE" = 'K_rankdir' ]; then

    mpiexec -n $NNP -of-proc $STDOUT ./${progdir}/${progbase} $CONF '' $ARGS
    res=$?
    if ((res != 0)); then
      echo "[Error] mpiexec -n $NNP -of-proc $STDOUT ./${progdir}/${progbase} $CONF '' $ARGS" >&2
      echo "        Exit code: $res" >&2
      exit $res
    fi

  else

    ( cd $progdir && mpiexec -n $NNP -of-proc $STDOUT ./$progbase $CONF '' $ARGS )
    res=$?
    if ((res != 0)); then 
      echo "[Error] mpiexec -n $NNP -of-proc $STDOUT ./$progbase $CONF '' $ARGS" >&2
      echo "        Exit code: $res" >&2
      exit $res
    fi

  fi

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

pdbash () {
#-------------------------------------------------------------------------------
# Submit bash parallel scripts according to nodefile
#
# Usage: pdbash NODEFILE PROC_OPT SCRIPT [ARGS]
#
#   NODEFILE  Name of nodefile (omit the directory $NODEFILE_DIR)
#   PROC_OPT  Options of using processes
#             all:  run the script in all processes listed in $NODEFILE
#             one:  run the script only in the first process and node in $NODEFILE
#   SCRIPT    Script (the working directory is set to $SCRP_DIR)
#   ARGS      Arguments passed into the program
#
# Other input variables:
#   $NODEFILE_DIR  Directory of nodefiles
#-------------------------------------------------------------------------------

if (($# < 2)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local NODEFILE="$1"; shift
local PROC_OPT="$1"; shift
local SCRIPT="$1"; shift
local ARGS="$@"

if [ -x "$TMPDAT/exec/pdbash" ]; then
  pdbash_exec="$TMPDAT/exec/pdbash"
elif [ -x "$COMMON_DIR/pdbash" ]; then
  pdbash_exec="$COMMON_DIR/pdbash"
else
  echo "[Error] $FUNCNAME: Cannot find 'pdbash' program." >&2
  exit 1
fi

if [ "$PROC_OPT" != 'all' ] && [ "$PROC_OPT" != 'one' ]; then
  echo "[Error] $FUNCNAME: \$PROC_OPT needs to be {all|one}." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

if [ "$MPI_TYPE" = 'sgimpt' ]; then

  if [ "$PROC_OPT" == 'all' ]; then
    local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
  elif [ "$PROC_OPT" == 'one' ]; then
    local HOSTLIST=$(head -n 1 ${NODEFILE_DIR}/${NODEFILE})
  fi
  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

  $MPIRUN -d $SCRP_DIR $HOSTLIST 1 $pdbash_exec $SCRIPT $ARGS
#  $MPIRUN -d $SCRP_DIR $HOSTLIST 1 bash $SCRIPT - $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $MPIRUN -d $SCRP_DIR $HOSTLIST 1 $pdbash_exec $SCRIPT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'openmpi' ]; then

  if [ "$PROC_OPT" == 'all' ]; then
    NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)
  elif [ "$PROC_OPT" == 'one' ]; then
    NNP=1
  fi

  $MPIRUN -np $NNP -hostfile ${NODEFILE_DIR}/${NODEFILE} -wdir $SCRP_DIR $pdbash_exec $SCRIPT $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $MPIRUN -np $NNP -hostfile ${NODEFILE_DIR}/${NODEFILE} -wdir $SCRP_DIR $pdbash_exec $SCRIPT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'K' ]; then

  if [ "$STG_TYPE" = 'K_rankdir' ]; then
    if [ "$PROC_OPT" == 'one' ]; then

      mpiexec -n 1 $pdbash_exec $SCRIPT $ARGS
      res=$?
      if ((res != 0)); then
        echo "[Error] mpiexec -n 1 $pdbash_exec $SCRIPT $ARGS" >&2
        echo "        Exit code: $res" >&2
        exit $res
      fi

    else

      mpiexec $pdbash_exec $SCRIPT $ARGS
      res=$?
      if ((res != 0)); then
        echo "[Error] mpiexec $pdbash_exec $SCRIPT $ARGS" >&2
        echo "        Exit code: $res" >&2
        exit $res
      fi

    fi
  else
    if [ "$PROC_OPT" == 'one' ]; then

      ( cd $SCRP_DIR && mpiexec -n 1 $pdbash_exec $SCRIPT $ARGS )
      res=$?
      if ((res != 0)); then
        echo "[Error] mpiexec -n 1 $pdbash_exec $SCRIPT $ARGS" >&2
        echo "        Exit code: $res" >&2
        exit $res
      fi

    else

      ( cd $SCRP_DIR && mpiexec $pdbash_exec $SCRIPT $ARGS )
      res=$?
      if ((res != 0)); then
        echo "[Error] mpiexec $pdbash_exec $SCRIPT $ARGS" >&2
        echo "        Exit code: $res" >&2
        exit $res
      fi

    fi
  fi

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

pdrun () {
#-------------------------------------------------------------------------------
# Return if it is the case to run parallel scripts, according to nodefile
#
# Usage: pdrun GROUP OPT
#
#   GROUP   Group of processes
#           all:     all processes
#           (group): process group #(group)
#   OPT     Options of the ways to pick up processes
#           all:  run the script in all processes in the group
#           alln: run the script in all nodes in the group, one process per node (default)
#           one:  run the script only in the first process in the group
#
# Other input variables:
#   MYRANK  The rank of the current process
#
# Exit code:
#   0: This process is used
#   1: This process is not used
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local GROUP="$1"; shift
local OPT="${1:-alln}"

#-------------------------------------------------------------------------------

local mynode=${proc2node[$((MYRANK+1))]}
if [ -z "$mynode" ]; then
  exit 1
fi

local res=1
local n

if [ "$GROUP" = 'all' ]; then

  if [ "$OPT" = 'all' ]; then
    exit 0
  elif [ "$OPT" = 'alln' ]; then
    res=0
    for n in $(seq $MYRANK); do
      if ((${proc2node[$n]} == mynode)); then
        res=1
        break
      fi
    done
  elif [ "$OPT" = 'one' ]; then
    if ((MYRANK == 0)); then
      exit 0
    fi
  fi

elif ((GROUP <= parallel_mems)); then

  local mygroup=${proc2group[$((MYRANK+1))]}
  local mygrprank=${proc2grpproc[$((MYRANK+1))]}

  if ((mygroup = GROUP)); then
    if [ "$OPT" = 'all' ]; then
      exit 0
    elif [ "$OPT" = 'alln' ]; then
      res=0
      for n in $(seq $((mygrprank-1))); do
        if ((${mem2node[$(((GROUP-1)*mem_np+n))]} == mynode)); then
          res=1
          break
        fi
      done
    elif [ "$OPT" = 'one' ]; then
      if ((${mem2node[$(((GROUP-1)*mem_np+1))]} == mynode)); then
        res=0
        for n in $(seq $((mygrprank-1))); do
          if ((${mem2node[$(((GROUP-1)*mem_np+n))]} == mynode)); then
            res=1
            break
          fi
        done
      fi
    fi
  fi

fi

exit $res

#-------------------------------------------------------------------------------
}

#===============================================================================

bdy_setting () {
#-------------------------------------------------------------------------------
# Calculate scale_init namelist settings for boundary files
#
# Usage: bdy_setting TIME FCSTLEN PARENT_LCYCLE [PARENT_FOUT] [PARENT_REF_TIME]
#
#   TIME
#   FCSTLEN
#   PARENT_LCYCLE
#   PARENT_FOUT
#   PARENT_REF_TIME
#
# Return variables:
#   $nbdy
#   $ntsteps
#   $ntsteps_skip
#   $ntsteps_total
#   $bdy_times[1...$nbdy]
#   $bdy_start_time
#   $parent_start_time
#
#  *Require source 'func_datetime' first.
#-------------------------------------------------------------------------------

if (($# < 4)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local TIME=$(datetime $1); shift
local FCSTLEN=$1; shift
local PARENT_LCYCLE=$1; shift
local PARENT_FOUT=${1:-$PARENT_LCYCLE}; shift
local PARENT_REF_TIME=${1:-$TIME}

if [ "$PARENT_FOUT" = '-' ]; then
  PARENT_FOUT=$FCSTLEN
fi
if [ "$PARENT_LCYCLE" = '-' ]; then
  PARENT_LCYCLE=$((FCSTLEN+PARENT_FOUT))
fi
PARENT_REF_TIME=$(datetime $PARENT_REF_TIME)

#-------------------------------------------------------------------------------
# compute $ntsteps

if ((PARENT_LCYCLE % PARENT_FOUT != 0)); then
  echo "[Error] $FUNCNAME: $PARENT_LCYCLE needs to be an exact multiple of $PARENT_FOUT." >&2
  exit 1
fi
ntsteps=$((PARENT_LCYCLE / PARENT_FOUT))

#-------------------------------------------------------------------------------
# compute $parent_start_time based on $PARENT_REF_TIME and $PARENT_LCYCLE

parent_start_time=$PARENT_REF_TIME
local parent_start_time_prev=$parent_start_time
while ((parent_start_time <= TIME)); do
  parent_start_time_prev=$parent_start_time
  parent_start_time=$(datetime $parent_start_time $PARENT_LCYCLE s)
done
parent_start_time=$parent_start_time_prev

while ((parent_start_time > TIME)); do
  parent_start_time=$(datetime $parent_start_time -${PARENT_LCYCLE} s)
done

#-------------------------------------------------------------------------------
# compute $bdy_start_time, $ntsteps_skip, and $ntsteps_total based on $parent_start_time and $PARENT_FOUT
# (assume $bdy_start_time <= $TIME)

ntsteps_skip=-1
bdy_start_time=$parent_start_time
while ((bdy_start_time <= TIME)); do
  bdy_start_time_prev=$bdy_start_time
  bdy_start_time=$(datetime $bdy_start_time $PARENT_FOUT s)
  ntsteps_skip=$((ntsteps_skip+1))
done
bdy_start_time=$bdy_start_time_prev

local ntsteps_total=$(((FCSTLEN-1)/PARENT_FOUT+2 + ntsteps_skip))

if ((bdy_start_time != TIME)); then
  if (($(datetime $bdy_start_time $(((ntsteps_total-1)*PARENT_FOUT)) s) < $(datetime $TIME $FCSTLEN s))); then
    ntsteps_total=$((ntsteps_total+1))
  fi
fi

#-------------------------------------------------------------------------------
# compute $nbdy and $bdy_times[1...$nbdy]

nbdy=1
bdy_times[1]=$parent_start_time
while ((ntsteps_total > ntsteps)); do
  nbdy=$((nbdy+1))
  bdy_times[$nbdy]=$(datetime ${bdy_times[$((nbdy-1))]} $PARENT_LCYCLE s)
  ntsteps_total=$((ntsteps_total-ntsteps))
done

if ((nbdy == 1)); then
  ntsteps=$ntsteps_total
fi

#echo "\$nbdy              = $nbdy" >&2
#echo "\$ntsteps           = $ntsteps" >&2
#echo "\$ntsteps_skip      = $ntsteps_skip" >&2
#echo "\$ntsteps_total     = $ntsteps_total" >&2
#echo "\$bdy_start_time    = $bdy_start_time" >&2
#echo "\$parent_start_time = $parent_start_time" >&2

#-------------------------------------------------------------------------------
}

#===============================================================================

job_submit_PJM () {
#-------------------------------------------------------------------------------
# Submit a PJM job.
#
# Usage: job_submit_PJM
#
#   JOBSCRP  Job script
#
# Return variables:
#   $jobid  Job ID monitered
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBSCRP="$1"

local rundir=$(dirname $JOBSCRP)
local scrpname=$(basename $JOBSCRP)

#-------------------------------------------------------------------------------

res=$(cd $rundir && pjsub $scrpname 2>&1)
echo $res

if [ -z "$(echo $res | grep '\[ERR.\]')" ]; then
  jobid=$(echo $res | grep 'submitted' | cut -d ' ' -f 6)
  if [ -z "$jobid" ]; then
    echo "[Error] $FUNCNAME: Error found when submitting a job." >&2
    exit 1
  fi
else
  echo "[Error] $FUNCNAME: Error found when submitting a job." >&2
  exit 1
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

job_end_check_PJM () {
#-------------------------------------------------------------------------------
# Check if a K-computer job has ended.
#
# Usage: job_end_check_PJM JOBID
#
#   JOBID  Job ID monitored
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBID="$1"

#-------------------------------------------------------------------------------

while true; do
  jobnum=$(pjstat $JOBID | sed -n '2p' | awk '{print $10}')
  if [[ "$jobnum" =~ ^[0-9]+$ ]]; then
    if ((jobnum == 0)); then
      break
    fi
  fi
  sleep 30s
done

#-------------------------------------------------------------------------------
}

#===============================================================================
