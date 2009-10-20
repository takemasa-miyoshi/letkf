#!/bin/sh
### Running on LINUX
#
# DIR
CDIR=`pwd`
cd ..
AFES=`pwd`
OBSDIR=$AFES/DATA/obs
DLDIR=$OBSDIR/downloads
SAVEDIR=$OBSDIR/prepbufr
GRABBUFR=/export/data/dvl0/miyoshi/downloads/grabbufr
BUFRLIB=/export/data/dvl0/miyoshi/downloads/BUFRLIB
F90=pgf90
FFLAGS="-byteswapio" #big-endian IO
# Initial time
IY=2008
IM=01
ID=01
IH=00
# Final time
EY=2008
EM=01
ED=01
EH=00
#
source $AFES/../common/timeinc.sh
#
# WKDIR
#
WKDIR=$OBSDIR/tmp
rm -rf $WKDIR
mkdir -p $WKDIR
cd $WKDIR
#
# mkdir
#
mkdir -p $DLDIR
mkdir -p $SAVEDIR
#
# build grabbufr
#
cp $GRABBUFR/grabbufr.f .
cp $GRABBUFR/spbufr.f .
$F90 -o grabbufr grabbufr.f spbufr.f
#
# build decoder
#
cp $AFES/../common/SFMT.f90 .
cp $AFES/../common/common.f90 .
cp $AFES/../common/common_obs.f90 .
cp $CDIR/dec_prepbufr.f90 .
$F90 $FFLAGS -o decoder SFMT.f90 common.f90 common_obs.f90 dec_prepbufr.f90 -L$BUFRLIB -lbufrlib
#
# Cycle run ### MAIN LOOP ###
#
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do
echo " >>>"
echo " >>> BEGIN COMPUTATION OF $IY/$IM/$ID/$IH"
echo " >>>"
#
# GET NCEP PREPBUFR
#
cd $DLDIR
DATAURL="http://nomads.ncdc.noaa.gov/data/gdas"
if test ! -f $DLDIR/$IY$IM$ID$IH.prepbufr.nr
then
wget $DATAURL/$IY$IM/$IY$IM$ID/gdas1.t${IH}z.prepbufr.nr
mv gdas1.t${IH}z.prepbufr.nr $IY$IM$ID$IH.prepbufr.nr
fi
#
# DECODE
#
cd $WKDIR
ln -s $DLDIR/$IY$IM$ID$IH.prepbufr.nr prepbufr.tmp
wc -c prepbufr.tmp | ./grabbufr prepbufr.tmp prepbufr.in
time ./decoder
#
# SAVE
#
mv fort.90 $SAVEDIR/obs$IY$IM$ID$IH.dat
#
# Date change ### MAIN LOOP END ###
#
TY=`timeinc6hr $IY $IM $ID $IH | cut -c1-4`
TM=`timeinc6hr $IY $IM $ID $IH | cut -c5-6`
TD=`timeinc6hr $IY $IM $ID $IH | cut -c7-8`
TH=`timeinc6hr $IY $IM $ID $IH | cut -c9-10`
IY=$TY
IM=$TY
ID=$TY
IH=$TY
done

