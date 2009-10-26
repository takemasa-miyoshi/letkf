#!/bin/sh
### Running on LINUX
#
# DATA SOURCE
#
DATASRC=2 #1:NOAA/NOMADS, 2:UCAR/DSS
#
# COMPILER
#
F90=pgf90
FFLAGS="-byteswapio" #big-endian IO
#
# TOOLS: grabbufr, BUFRLIB
#
GRABBUFR=/export/data/dvl0/miyoshi/downloads/grabbufr
BUFRLIB=/export/data/dvl0/miyoshi/downloads/BUFRLIB
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
# DIR
#
CDIR=`pwd`
cd ..
AFES=`pwd`
OBSDIR=$AFES/DATA/obs
DLDIR=$OBSDIR/downloads
SAVEDIR=$OBSDIR/prepbufr
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
# Settings for DATASRC
#
if test $DATASRC -eq 1
then
DATAURL="http://nomads.ncdc.noaa.gov/data/gdas"
WGET="wget"
DLDIR=$OBSDIR/downloads/nomads
SAVEDIR=$OBSDIR/prepbufr/nomads
mkdir -p $DLDIR
mkdir -p $SAVEDIR
else
EMAIL=$1
PASSWD=$2
if test $# -ne 2
then
echo "USAGE: $0 EMAIL PASSWD"
exit
fi
DATAURL="https://dss.ucar.edu/datazone/dsszone/ds337.0"
WGET="wget --no-check-certificate"
DLDIR=$OBSDIR/downloads/ucar
SAVEDIR=$OBSDIR/prepbufr/ucar
mkdir -p $DLDIR
mkdir -p $SAVEDIR
cd $DLDIR
rm -f auth.dss_ucar_edu
OPT1='-O /dev/null --save-cookies auth.dss_ucar_edu --post-data'
OPT2="email=${EMAIL}&passwd=${PASSWD}&action=login"
$WGET $OPT1="$OPT2" https://dss.ucar.edu/cgi-bin/login
OPT="-N --load-cookies auth.dss_ucar_edu"
cd $WKDIR
fi
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
if test $DATASRC -eq 1
then
DATAF=$IY$IM$ID$IH.prepbufr.nr
else
DATAF=prepbufr.gdas.$IY$IM$ID$IH.wo40
fi
if test ! -f $DATAF
then
if test $DATASRC -eq 1
then
$WGET $DATAURL/$IY$IM/$IY$IM$ID/gdas1.t${IH}z.prepbufr.nr
mv gdas1.t${IH}z.prepbufr.nr $IY$IM$ID$IH.prepbufr.nr
else
$WGET $OPT $DATAURL/prepbufr.$IY$IM$ID.wo40.tar.gz
tar zxvf prepbufr.$IY$IM$ID.wo40.tar.gz
rm prepbufr.$IY$IM$ID.wo40.tar.gz
mv $IY$IM$ID.wo40/* .
rmdir $IY$IM$ID.wo40
rm -rf $IY$IM$ID.wo40.tar.gz
fi
fi
#
# DECODE
#
cd $WKDIR
ln -s $DLDIR/$DATAF prepbufr.tmp
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
IM=$TM
ID=$TD
IH=$TH
done
rm -f $DLDIR/auth.dss_ucar_edu
echo "NORMAL END"
