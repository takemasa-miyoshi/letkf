#!/bin/sh
#=======================================================================
# ensfcst.sh
#   This script runs the SPEEDY model with subdirectory $NODE
#=======================================================================
set -e
### input for this shell
SPEEDY=$1
OUTPUT=$2
YMDH=$3
TYMDH=$4
MEM=$5
NODE=$6
###
if test 5$6 -eq 5
then
echo "ERROR in ensfcst.sh"
exit
fi
### run
rm -rf $NODE
mkdir $NODE
cp imp.exe $NODE
SB=$SPEEDY/model/data/bc/t30/clim
SC=$SPEEDY/model/data/bc/t30/anom
ln -s $SB/orog_lsm_alb.t30.grd         $NODE/fort.20
ln -s $SB/sst_8190clim.t30.sea.grd     $NODE/fort.21
ln -s $SB/seaice_8190clim.t30.sea.grd  $NODE/fort.22
ln -s $SB/skt_8190clim.t30.land.grd    $NODE/fort.23
ln -s $SB/sndep_8190clim.t30.land.grd  $NODE/fort.24
ln -s $SB/veget.t30.land.grd           $NODE/fort.25
ln -s $SB/soilw_8190clim.t30.land.grd  $NODE/fort.26
cp    $SC/sst_anom_7990.t30.grd        $NODE/fort.30
### run
cd $NODE
ln -fs $OUTPUT/anal/$MEM/$YMDH.grd fort.90
FORT2=1111
echo $FORT2 > fort.2
echo $YMDH | cut -c1-4 >> fort.2
echo $YMDH | cut -c5-6 >> fort.2
echo $YMDH | cut -c7-8 >> fort.2
echo $YMDH | cut -c9-10 >> fort.2
./imp.exe > out.lis 2> out.lis.2
### output
mv ${YMDH}.grd $OUTPUT/anal_f/$MEM
mv ${YMDH}_p.grd $OUTPUT/anal_f/$MEM
mv ${TYMDH}.grd $OUTPUT/gues/$MEM
mv ${TYMDH}_p.grd $OUTPUT/gues/$MEM
exit 0
