#!/bin/bash
#===============================================================================
#
#  Main settings for SCALE-LETKF scripts
#
#===============================================================================

PRESET='K_rankdir'                          # 'K' or 'K_rankdir' or 'K_micro'

#===============================================================================

DIR="$(cd "$(pwd)/.." && pwd)"                  # Root directory of the SCALE-LETKF

DDIR="$(cd "$(pwd)/../../../.." && pwd)"    # Root directory of the test suite

OUTDIR="$DDIR/exp/testcase_45km_4p_l72"     # Directory for SCALE-LETKF output

#===============================================================================
# Location of model/data files

MODELDIR="$DIR/../../bin/${SCALE_SYS}"                               # Directory of the SCALE model executables
DATADIR="$DDIR/database"                                             # Directory of the SCALE database

DATA_TOPO=                                                           # Directory of the prepared topo files
DATA_TOPO_BDY_SCALE=
DATA_LANDUSE=                                                        # Directory of the prepared landuse files
DATA_BDY_SCALE=                                                      # Directory of the boundary data in SCALE history format (parent domain)
DATA_BDY_SCALE_PREP=                                                 # Directory of the prepared SCALE boundary files
DATA_BDY_WRF="$DDIR/ncepfnl/wrfout"                                  # Directory of the boundary data in WRF format
DATA_BDY_NICAM=                                                      # Directory of the boundary data in NICAM format (not finished)

OBS="$DDIR/obs/prepbufr_obs_eastasia"                                # Observation data in LETKF observation format
OBSNCEP=

#===============================================================================
# model/data file options

TOPO_FORMAT='GTOPO30'   # 'prep': Use prepared topo files in $DATA_TOPO
                        # 'GTOPO30' (requires compatible 'config.nml.scale_pp')
                        # 'DEM50M'  (requires compatible 'config.nml.scale_pp')

LANDUSE_FORMAT='GLCCv2' # 'prep': Use prepared landuse files in $DATA_LANDUSE
                        # 'GLCCv2' (requires compatible 'config.nml.scale_pp')
                        # 'LU100M' (requires compatible 'config.nml.scale_pp')
LANDUSE_UPDATE=0        # 0: Time-invariant landuse files
                        # 1: Time-variant landuse files

BDY_FORMAT=2            # 0: SCALE boundary files (with exactly same domain settings; do not need additional preprocessing)
                        # 1: SCALE history (requires compatible 'config.nml.scale_init')
                        # 2: WRF           (requires compatible 'config.nml.scale_init')
                        # 3: NICAM         (requires compatible 'config.nml.scale_init')
BDY_SCALE_DIR='hist'    # Directory name of the SCALE history files when $BDY_FORMAT = 1

BDY_ENS=0               # 0: Fixed boundary files for all memebers
                        # 1: Ensemble boundary files
BDY_ROTATING=0          # 0: Use a same series of boundary files for all initial time
                        # 1: Use different series of boundary files for different initial time

BDYINT=21600
BDYCYCLE_INT=

PARENT_REF_TIME=

OCEAN_INPUT=0           # 0: No ocean input (use cycling ocean variables)
                        # 1: Update the ocean variables every cycle

OCEAN_FORMAT=99         # 0: SCALE init files (with exactly same domain settings; do not need additional preprocessing)
                        # 99: From the same file as used in generating the boundary conditions ($BDY_FORMAT)

OBSNUM=1
OBSNAME[1]=obs                                    
OBSNAME[2]=radar

#===============================================================================
# Cycling settings

WINDOW_S=10800     # SCALE forecast time when assimilation window starts (second)
WINDOW_E=32400     # SCALE forecast time when assimilation window ends (second)
LCYCLE=21600       # Length of a GFS-LETKF cycle (second)
LTIMESLOT=3600     # Timeslot interval for 4D-LETKF (second)

#===============================================================================
# Parallelization settings

MEMBER=10          # Ensemble size

NNODES=44          # Number of nodes
PPN=1              # Number of processes per node

THREADS=8          # Number of threads per process

SCALE_NP=4         # Number of processes to run SCALE

BGJOB_INT='0.1s'   # Interval of multiple background job submissions

#ENABLE_SET=1       ######

#===============================================================================
# Temporary directories to store runtime files

DATA_BDY_TMPLOC=1           # Location of the temporary directory for DATA_BDY
                            #  1: in $TMPDAT
                            #  2: in $TMPOUT

DISK_MODE_TOPO_LANDUSE_DB=2

DISK_MODE_DATA_BDY=2

DISK_MODE_DATA_TOPO=2
DISK_MODE_DATA_LANDUSE=2
DISK_MODE_DATA_BDYPREP=2

ONLINE_STGOUT=0             # Stage out right after each cycle (do not wait until the end of the job)?
                            #  0: No
                            #  1: Yes

#SIMPLE_STGOUT=1             # Stage out the entire 'out' temporary directory instead of specifying each file?
#                            #  0: No
#                            #  1: Yes

SYSNAME="$(basename $OUTDIR)"                # A unique name in the machine
TMPSUBDIR="scale-letkf_${SYSNAME}" # (used to identify multiple runs in the same time)

TMP="/scratch/$(id -ng)/gylien/$TMPSUBDIR" # Temporary directory shared among all nodes
TMPS="$DIR/tmp/$TMPSUBDIR"  # Temporary directory only on the server node
#TMPL=

CLEAR_TMP=0                 # Clear temporary directories after the completion of job?
                            #  0: No
                            #  1: Yes

#===============================================================================
# Environmental settings

MPIBIN=$(dirname $(which mpiexec))
MPIRUN="$MPIBIN/mpiexec"

SCP='cp -L'
SCP_HOSTPREFIX=''
#SCP="scp -q"
#SCP_HOSTPREFIX="XXXX:"

SCP_THREAD=8
TAR_THREAD=8

PYTHON="python"

#BUFRBIN=

#===============================================================================
# Machine-independent source file

. config.rc

#===============================================================================
