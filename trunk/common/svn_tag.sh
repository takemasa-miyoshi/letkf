#!/bin/sh
REPO=https://miyoshi.googlecode.com/svn
TAG=alera2.0
USER=takemiyo
COMMONSRC="SFMT.f90 common.f90 common_mpi.f90 common_obs.f90 common_mtx.f90 common_letkf.f90 netlib.f"
for SRC in $COMMONSRC
do
svn copy $REPO/trunk/common/$SRC $REPO/tags/$TAG/common/$SRC -m "tagging $TAG" --username $USER
done
svn copy $REPO/trunk/afes $REPO/tags/$TAG/afes -m "tagging $TAG" --username $USER

