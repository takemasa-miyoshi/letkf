#!/bin/sh
#  plot bred vectors
set -e
TIME=0020
MEM=mean
VAR=v
XMIN=-770
XMAX=0
Z=40
CMAX=1
CMIN=-2
CINT=0.1
CDIR=`pwd`
cd ../..
ROMS=`pwd`
NATURE=$ROMS/DATA/nature
EXP=letkf001
Z1=`expr $Z - 1` || test 1 -eq 1
cd $CDIR
cat << EOF > tmp.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;wks=gsn_open_wks("x11","gsun02n")
wks=gsn_open_wks("eps","out")
;gsn_define_colormap(wks,"testcmap")
gsn_define_colormap(wks,"rainbow")
;gsn_define_colormap(wks,"gui_default")
;gsn_define_colormap(wks,"gsltod")
plot=new(3,graphic)
data1=addfile("$NATURE/$TIME.nc","r")
data2=addfile("$ROMS/DATA/$EXP/gues/$MEM/$TIME.nc","r")
data3=addfile("$ROMS/DATA/$EXP/anal/$MEM/$TIME.nc","r")
xlab=fspan(-770,21,114)
;ylab=fspan(0,791,114)
ylab=fspan(0,784,113)

z1=data1->$VAR(0,$Z1,:,:)
;z1=data1->zeta(0,:,:)
z1&eta_v=ylab
z1&xi_rho=xlab
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
;z2=z2-z1
;z2=data2->zeta(0,:,:)
;printVarSummary(z2)
z2&eta_v=ylab
z2&xi_rho=xlab
res@tiMainString="GUES(${MEM})"
;res@cnMaxLevelValF=1
;res@cnMinLevelValF=0
;res@cnLevelSpacingF=0.1
plot(1)=gsn_csm_contour(wks,z2,res)

zz=data3->$VAR(0,$Z1,:,:)
z3=zz
;z3=z3-z1
;z3=z3*10
;z3=data3->zeta(0,:,:)
z3&eta_v=ylab
z3&xi_rho=xlab
;printVarSummary(z3)
res@tiMainString="ANAL(${MEM})"
;res@cnMaxLevelValF=1
;res@cnMinLevelValF=-1
;res@cnLevelSpacingF=0.1
plot(2)=gsn_csm_contour(wks,z3,res)

gsn_panel(wks,plot,(/3,1/),False)
end
EOF
ncl tmp.ncl
rm tmp.ncl
display out.eps
echo "NORMAL END"
