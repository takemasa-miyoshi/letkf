#!/bin/sh
#=======================================================================
# ensfcst.sh
#   This script runs the ROMS model with subdirectory $NODE
#=======================================================================
set -e
### input for this shell
ROMS=$1
OUTPUT=$2
YMDH=$3
TYMDH=$4
MEM=$5
NODE=$6
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
ln -fs $OUTPUT/anal/$MEM/${YMDH}_rst.nc input.nc
./roms > roms.log
rm input.nc
mv output.0001.nc grd.nc
./chtimestep
mv grd.nc $OUTPUT/gues/$MEM/${TYMDH}_rst.nc
exit 0
