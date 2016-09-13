#!/bin/bash
#===============================================================================
#
#  Wrap fcst.sh in a torque job script and run it.
#
#  November 2014, created,                 Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_torque.sh [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP TIME_LIMIT]
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

#. src/func_distribute.sh
. src/func_datetime.sh
#. src/func_util.sh
. src/func_$myname1.sh

#-------------------------------------------------------------------------------

setting "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}"

jobscrp="${myname1}_job.sh"

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $(basename $0) $@"
echo

for vname in DIR OUTDIR DATA_TOPO DATA_LANDUSE DATA_BDY DATA_BDY_WRF OBS OBSNCEP MEMBER NNODES PPN \
             FCSTLEN FCSTOUT EFSOFLEN EFSOFOUT OUT_OPT \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP; do
  printf '  %-10s = %s\n' $vname "${!vname}"
done

echo

#===============================================================================
# Creat a job script

echo "[$(datetime_now)] Create a job script '$jobscrp'"

cat > $jobscrp << EOF
#!/bin/sh
#PBS -N ${myname1}_${SYSNAME}
#PBS -l "nodes=${NNODES}:ppn=${PPN}"
##PBS -l "select=${NNODES}:ncpus=$THREADS"
#PBS -l "walltime=${TIME_LIMIT}"
##PBS -W umask=027

ulimit -s unlimited

HOSTLIST=\$(cat \$PBS_NODEFILE | sort | uniq)
HOSTLIST=\$(echo \$HOSTLIST | sed 's/  */,/g')
export MPI_UNIVERSE="\$HOSTLIST $((PPN*THREADS))"

export OMP_NUM_THREADS=${THREADS}
#export PARALLEL=${THREADS}

export FORT_FMT_RECL=400

cd \$PBS_O_WORKDIR

rm -f machinefile
cp -f \$PBS_NODEFILE machinefile

./${myname1}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP"
EOF

echo "[$(datetime_now)] Run ${myname1} job on PBS"
echo



exit



job_submit_PJM $jobscrp
echo

job_end_check_PJM $jobid

#===============================================================================

echo
echo "[$(datetime_now)] Finish $(basename $0) $@"

exit 0
