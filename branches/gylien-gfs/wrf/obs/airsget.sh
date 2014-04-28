#/bin/sh
#=======================================================================
# [PURPOSE:]
#   download AIRS IR, AMSU, HSB retrievals -> LETKF format output
#
# [REQUIRED LIBRARIES:]
# 1. HDF-related: hdfeos, hdf5, hdf
# 2. Other standard library (usually in /usr/lib): libjpeg.a libsz.a, libz.a
#
# [PREPARING AIRS DATA URL:]
# 1. Go to http://disc.sci.gsfc.nasa.gov/AIRS/data-holdings/by-data-product
# 2. Go down to AIRS Level-2 Products, and find a row of "AIRX2RET"
# 3. Click a link of "Search" on the rightmost column "GES DISC Data Access"
# 4. Input Time Span and Location (optional) and click "Search GES-DISK"
# 5. "Add Selected Files To Cart"
# 6. "Continue to Shopping Cart" and "Checkout". Note: this is free.
# 7. Click "URL List (Data)"
# 8. Save the file list to "airslist.txt", and remove the final line (README)
#
# [AUTHOR:] T. Miyoshi  March 3, 2011, College Park, MD, USA
#=======================================================================
set -e
F90=pgf90
FFLAGS="-byteswapio"
HDFEOS=/usr/local/hdfeos/lib/linux
HDF5=/usr/local/hdf5/lib
HDF=/usr/local/hdf64/lib
#
# DIR
#
CDIR=`pwd`
cd ..
WRF=`pwd`
OBSDIR=$WRF/DATA/obs
DLDIR=$OBSDIR/downloads
SAVEDIR=$OBSDIR/airsret
mkdir -p $DLDIR/airs
mkdir -p $SAVEDIR
EXIST=`ls $SAVEDIR`
if test -n "$EXIST"
then
echo " >> In order to avoid duplication of sounding data, $SAVEDIR is going to be cleaned up. Is it OK? (y/n)"
read ANS
if test $ANS != 'y'
then
echo " >> Please check $SAVEDIR"
exit
else
echo " >> Clearning $SAVEDIR"
rm -rf $SAVEDIR/*
fi
fi
#
# WKDIR
#
WKDIR=$OBSDIR/tmp
rm -rf $WKDIR
mkdir -p $WKDIR
cd $WKDIR
#
# Make decoder
#
echo " >> Making decoder.."
cp $WRF/../common/SFMT.f90 .
cp $WRF/../common/common.f90 .
cp $WRF/common/common_wrf.f90 .
cp $WRF/common/common_obs_wrf.f90 .
cp $CDIR/dec_airsret.f90 .
$F90 $FFLAGS -R$HDF5 -o decoder SFMT.f90 common.f90 common_wrf.f90 common_obs_wrf.f90 dec_airsret.f90 -L$HDFEOS -lhdfeos -lGctp -L$HDF5 -lhdf5 -L$HDF -lmfhdf -ldf -lz -ljpeg -lsz
#
# Download and decode
#
cp $CDIR/airslist.txt .
NFILES=`wc -l airslist.txt | awk '{print $1}'`
N=1
while test $N -le $NFILES
do
URL=`cat airslist.txt | head -n $N | tail -n 1`
FILE=`basename $URL`
if test -f $DLDIR/airs/$FILE
then
echo " >> $FILE exists. skip downloading."
cp $DLDIR/airs/$FILE .
else
echo " >> downloading.. $FILE"
wget $URL
cp $FILE $DLDIR/airs/$FILE
fi
echo $FILE | ./decoder
DFILES=`ls *.dat`
for DFILE in $DFILES
do
cat $DFILE >> $SAVEDIR/$DFILE
rm $DFILE
done
N=`expr $N + 1`
done

echo "NORMAL END"
