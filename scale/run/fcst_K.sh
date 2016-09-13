#!/bin/bash
#===============================================================================
#
#  Wrap fcst.sh in a K-computer job script and run it.
#
#  October 2014, created,                 Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_K.sh [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP TIME_LIMIT]
#
#===============================================================================

cd "$(dirname "$0")"
myname1='fcst'

#===============================================================================
# Configuration

. config.main
res=$? && ((res != 0)) && exit $res
. config.$myname1
res=$? && ((res != 0)) && exit $res

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_$myname1.sh

#-------------------------------------------------------------------------------

if ((TMPDAT_MODE == 1 || TMPRUN_MODE == 1 || TMPOUT_MODE == 1)); then
  echo "[Error] $0: When using a regular resource group," >&2
  echo "        \$TMPDAT_MODE, \$TMPRUN_MODE, \$TMPOUT_MODE all need to be 2 or 3." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $(basename $0) $@"

setting "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}"

echo
print_setting
echo

#===============================================================================
# Create and clean the temporary directory

echo "[$(datetime_now)] Create and clean the temporary directory"
 
safe_init_tmpdir $TMPS

#===============================================================================
# Determine the distibution schemes

echo "[$(datetime_now)] Determine the distibution schemes"

# K computer
NNODES_real=$NNODES
PPN_real=$PPN
NNODES=$((NNODES*PPN))
PPN=1

declare -a node
declare -a node_m
declare -a name_m
declare -a mem2node
declare -a mem2proc
declare -a proc2node
declare -a proc2group
declare -a proc2grpproc

safe_init_tmpdir $TMPS/node
distribute_fcst "$MEMBERS" $CYCLE - $TMPS/node

if ((CYCLE == 0)); then
  CYCLE=$cycle_auto
fi

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

STAGING_DIR="$TMPS/staging"

safe_init_tmpdir $STAGING_DIR
staging_list

#-------------------------------------------------------------------------------

cp $SCRP_DIR/config.main $TMPS

echo "SCRP_DIR=\"\$(pwd)\"" >> $TMPS/config.main
echo "NODEFILE_DIR=\"\$(pwd)/node\"" >> $TMPS/config.main

echo "NNODES=$NNODES" >> $TMPS/config.main
echo "PPN=$PPN" >> $TMPS/config.main
echo "NNODES_real=$NNODES_real" >> $TMPS/config.main
echo "PPN_real=$PPN_real" >> $TMPS/config.main

echo "PARENT_REF_TIME=$PARENT_REF_TIME" >> $TMPS/config.main

echo "RUN_LEVEL='K'" >> $TMPS/config.main

#===============================================================================
# Creat a job script

jobscrp="${myname1}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

if ((NNODES_real > 36864)); then
  rscgrp="huge"
elif ((NNODES_real > 384)); then
  rscgrp="large"
else
  rscgrp="small"
fi

cat > $jobscrp << EOF
#!/bin/sh
#PJM -N ${myname1}_${SYSNAME}
#PJM -s
#PJM --rsc-list "node=${NNODES_real}"
#PJM --rsc-list "elapse=${TIME_LIMIT}"
#PJM --rsc-list "rscgrp=${rscgrp}"
##PJM --rsc-list "node-quota=29G"
##PJM --mpi "shape=${NNODES_real}"
#PJM --mpi "proc=$NNODES"
#PJM --mpi assign-online-node
#PJM --stg-transfiles all
EOF

if [ "$STG_TYPE" = 'K_rankdir' ]; then
  echo "#PJM --mpi \"use-rankdir\"" >> $jobscrp
fi

bash $SCRP_DIR/src/stage_K.sh $STAGING_DIR $myname1 >> $jobscrp

cat >> $jobscrp << EOF

. /work/system/Env_base_1.2.0-20-1
export OMP_NUM_THREADS=${THREADS}
export PARALLEL=${THREADS}

./${myname1}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" || exit \$?
EOF

#===============================================================================
# Check the staging list

echo "[$(datetime_now)] Run pjstgchk"
echo

pjstgchk $jobscrp
res=$? && ((res != 0)) && exit $res
echo

#===============================================================================
# Run the job

echo "[$(datetime_now)] Run ${myname1} job on PJM"
echo

job_submit_PJM $jobscrp
echo

job_end_check_PJM $jobid

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

n=0
nmax=12
while [ ! -s "${myname1}_${SYSNAME}.i${jobid}" ] && ((n < nmax)); do
  n=$((n+1))
  sleep 5s
done

mkdir -p $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/config.main $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/config.${myname1} $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/config.nml.* $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/${myname1}_job.sh $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f ${myname1}_${SYSNAME}.o${jobid} $OUTDIR/exp/${jobid}_${myname1}_${STIME}/job.o
cp -f ${myname1}_${SYSNAME}.e${jobid} $OUTDIR/exp/${jobid}_${myname1}_${STIME}/job.e
cp -f ${myname1}_${SYSNAME}.i${jobid} $OUTDIR/exp/${jobid}_${myname1}_${STIME}/job.i
cp -f ${myname1}_${SYSNAME}.s${jobid} $OUTDIR/exp/${jobid}_${myname1}_${STIME}/job.s
( cd $SCRP_DIR ; git log -1 --format="SCALE-LETKF version %h (%ai)" > $OUTDIR/exp/${jobid}_${myname1}_${STIME}/version )
( cd $MODELDIR ; git log -1 --format="SCALE       version %h (%ai)" >> $OUTDIR/exp/${jobid}_${myname1}_${STIME}/version )

finalization

if ((CLEAR_TMP == 1)); then
  safe_rm_tmpdir $TMPS
fi

#===============================================================================

echo "[$(datetime_now)] Finish $(basename $0) $@"

exit 0
