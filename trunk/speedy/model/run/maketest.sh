#!/bin/sh
# maketest
#   For compile only
set -e

cd ../..
SPEEDY=`pwd`
#
# Work directory
#
TMPDIR=$SPEEDY/model/tmp
rm -rf $TMPDIR
mkdir -p $TMPDIR
cd $TMPDIR
#
# BUILD SPEEDY MODEL
#
echo '>>>BEGIN BUILDING SPEEDY MODEL'
cp $SPEEDY/model/source/makefile ./
cp $SPEEDY/model/source/*.h ./
cp $SPEEDY/model/source/*.f ./
cp $SPEEDY/model/source/*.s ./

mv par_horres_t30.h atparam.h
mv par_verres.h atparam1.h

cp $SPEEDY/model/ver32.input/cls_*.h ./
cp $SPEEDY/model/ver32.input/inpfiles.s ./

cp $SPEEDY/model/update/*.h ./
cp $SPEEDY/model/update/*.f ./
cp $SPEEDY/model/update/makefile ./

make imp.exe
echo '>>>END BUILDING SPEEDY MODEL'

