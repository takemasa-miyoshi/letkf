#!/bin/sh
#  plot bred vectors
set -e
OBS=5x20
EXP=LETKF20H80V20
TDV=TDVAR
TIME=2004013118
MEM=mean
VARLIST="u v temp salt"
ZLIST="40 30 20 10"
OUTPUT=figs/$OBS/${EXP}-$TDV
mkdir -p $OUTPUT
CDIR=`pwd`
cd ../..
ROMS=`pwd`
NATURE=$ROMS/DATA/nature
cd $CDIR
for VAR in $VARLIST
do
case $VAR in
  "u"    ) XMIN=-763;XLAB="fspan(-763,21,113)";YLAB="fspan(0,791,114)";XI="xi_u";ETA="eta_rho";CMAX=0.4;CMIN=-0.4;CINT=0.05;DMAX=0.5;DINT=0.05;;
  "v"    ) XMIN=-770;XLAB="fspan(-770,21,114)";YLAB="fspan(0,784,113)";XI="xi_rho";ETA="eta_v";CMAX=0.4;CMIN=-0.4;CINT=0.05;DMAX=0.5;DINT=0.05;;
  "temp" ) XMIN=-770;XLAB="fspan(-770,21,114)";YLAB="fspan(0,791,114)";XI="xi_rho";ETA="eta_rho";CMAX=20.0;CMIN=15.0;CINT=0.5;DMAX=1.0;DINT=0.1;;
  "salt" ) XMIN=-770;XLAB="fspan(-770,21,114)";YLAB="fspan(0,791,114)";XI="xi_rho";ETA="eta_rho";CMAX=34.0;CMIN=33.0;CINT=0.1;DMAX=0.5;DINT=0.05;;
  *      ) echo "ABEND";exit;;
esac
XMAX=0
for Z in $ZLIST
do
Z1=`expr $Z - 1` || test 1 -eq 1
case $VAR$Z in
  "u10"    ) CMAX=0.1;CMIN=-0.1;CINT=0.01;DMAX=0.1;DINT=0.01;;
  "v10"    ) CMAX=0.1;CMIN=-0.1;CINT=0.01;DMAX=0.1;DINT=0.01;;
  "u20"    ) CMAX=0.1;CMIN=-0.1;CINT=0.01;DMAX=0.1;DINT=0.01;;
  "v20"    ) CMAX=0.1;CMIN=-0.1;CINT=0.01;DMAX=0.1;DINT=0.01;;
  "temp20" ) CMAX=11.0;CMIN=7.0;CINT=0.1;;
  "temp10" ) CMAX=5.0;CMIN=4.0;CINT=0.1;DMAX=0.2;DINT=0.02;;
  "salt10" ) CMAX=34.4;CMIN=34.3;CINT=0.01;DMAX=0.03;DINT=0.003;;
esac
cat << EOF > tmp.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;wks=gsn_open_wks("x11","gsun02n")
wks=gsn_open_wks("eps","out")
;gsn_define_colormap(wks,"testcmap")
gsn_define_colormap(wks,"rainbow")
;gsn_define_colormap(wks,"gui_default")
;gsn_define_colormap(wks,"gsltod")
plot=new(4,graphic)
data1=addfile("$NATURE/${TIME}_rst.nc","r")
data2=addfile("$ROMS/DATA/$OBS/$EXP/anal/$MEM/${TIME}_rst.nc","r")
data3=addfile("$ROMS/DATA/$OBS/$TDV/anal/${TIME}_da.nc","r")
data4=addfile("$ROMS/DATA/$OBS/$EXP/infl_mul/${TIME}_rst.nc","r")
xlab=$XLAB
ylab=$YLAB

z1=data1->$VAR(0,$Z1,:,:)
;z1=data1->zeta(0,:,:)
z1&$ETA=ylab
z1&$XI=xlab
;printVarSummary(z1)
res=True
res@tiMainString="NATURE"
res@gsnDraw=False
res@gsnFrame=False
res@gsnSpreadColors=True
res@cnFillOn=True
res@cnLinesOn=False
res@cnLevelSelectionMode="ManualLevels"
res@cnMaxLevelValF=$CMAX
res@cnMinLevelValF=$CMIN
res@cnLevelSpacingF=$CINT
res@cnLineLabelsOn=False
res@lbLabelAutoStride=True
res@lbOrientation="Vertical"
res@trXMinF=$XMIN
res@trXMaxF=$XMAX
plot(0)=gsn_csm_contour(wks,z1,res)

z=data2->$VAR(0,$Z1,:,:)
z2=z
delete(z)
;z2=z2-z1
;z2=data2->zeta(0,:,:)
;printVarSummary(z2)
z2&$ETA=ylab
z2&$XI=xlab
res@tiMainString="$EXP"
;res@cnMaxLevelValF=1
;res@cnMinLevelValF=0
;res@cnLevelSpacingF=0.1
plot(1)=gsn_csm_contour(wks,z2,res)

z=data3->$VAR(0,$Z1,:,:)
z3=dble2flt(z)
delete(z)
;z3=z3-z1
;z3=z3*10
;z3=data3->zeta(0,:,:)
z3&$ETA=ylab
z3&$XI=xlab
;printVarSummary(z3)
res@tiMainString="3DVAR"
;res@cnMaxLevelValF=1
;res@cnMinLevelValF=-1
;res@cnLevelSpacingF=0.1
plot(2)=gsn_csm_contour(wks,z3,res)

z=data1->$VAR(0,$Z1,:,:)
z4=z
z4=z2-z3
;delete(z)
z4&$ETA=ylab
z4&$XI=xlab
;printVarSummary(z4)
res@tiMainString="DIFF ${EXP}-3DVAR"
res@cnMaxLevelValF=$DMAX
res@cnMinLevelValF=-$DMAX
res@cnLevelSpacingF=$DINT
plot(3)=gsn_csm_contour(wks,z4,res)

gsn_panel(wks,plot,(/4,1/),False)
end
EOF
ncl tmp.ncl
rm tmp.ncl
mv out.eps $OUTPUT/plt_${VAR}_$Z.eps
done
done
echo "NORMAL END"
