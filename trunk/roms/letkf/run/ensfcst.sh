#!/bin/sh
#=======================================================================
# ensfcst.sh
#   This script runs the ROMS model with subdirectory $NODE
#=======================================================================
set -e
### input for this shell
ROMS=$1
OUTPUT=$2
TIME=$3
MEM=$4
NODE=$5
###
if test 5$5 -eq 5
then
echo "ERROR in ensfcst.sh"
exit
fi
### run
rm -rf $NODE
mkdir $NODE
cp $ROMS/model/src/roms $NODE
cp $ROMS/model/bc/*.nc $NODE
cp $ROMS/ncio/chtimestep $NODE
cp $ROMS/model/updates/roms.in.6hr $NODE/roms.in
### run
cd $NODE
ln -fs $OUTPUT/anal/$MEM/$TIME.nc input.nc
./roms > roms.log
rm input.nc
mv output.0001.nc grd.nc
./chtimestep
TIME=`expr $TIME + 1`
if test $TIME -lt 1000
then
TIME=0$TIME
fi
if test $TIME -lt 100
then
TIME=0$TIME
fi
if test $TIME -lt 10
then
TIME=0$TIME
fi
mv grd.nc $OUTPUT/gues/$MEM/$TIME.nc
exit 0
