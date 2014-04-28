#!/bin/sh
#  plot bred vectors
set -e
TIME=0221
#TIME=0701
MEM=001
XMIN=-770
XMAX=0
Z=5
CDIR=`pwd`
cd ../..
ROMS=`pwd`
NATURE=$ROMS/DATA/nature
BVEXP=bv_t5_5E-3
Z1=`expr $Z - 1` || test 1 -eq 1
#BVEXP=breeding
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
data1=addfile("$NATURE/$TIME.nc","r")
data2=addfile("$ROMS/DATA/$BVEXP/bv/$MEM/$TIME.nc","r")
xlab=fspan(-770,21,114)
ylab=fspan(0,791,114)

z1=data1->temp(0,$Z1,:,:)
;z1=data1->zeta(0,:,:)
z1&eta_rho=ylab
z1&xi_rho=xlab
;printVarSummary(z1)
res1=True
res1@tiMainString="NATURE(shade), BV(contour)"
res1@gsnDraw=False
res1@gsnFrame=False
res1@gsnSpreadColors=True
res1@cnFillOn=True
res1@cnLinesOn=False
res1@cnLevelSelectionMode="ManualLevels"
res1@cnMaxLevelValF=2.45
res1@cnMinLevelValF=2.15
res1@cnLevelSpacingF=0.005
;res1@cnMaxLevelValF=20
;res1@cnMinLevelValF=10
;res1@cnLevelSpacingF=0.1
res1@cnLineLabelsOn=False
res1@lbLabelAutoStride=True
res1@lbOrientation="Vertical"
res1@trXMinF=$XMIN
res1@trXMaxF=$XMAX
plot1=gsn_csm_contour(wks,z1,res1)

z2=data2->temp(0,$Z1,:,:)
;z2=data2->zeta(0,:,:)
z2=z2/max(abs(z2))
z2&eta_rho=ylab
z2&xi_rho=xlab
;printVarSummary(z2)
res2=True
res2@tiMainString=""
res2@gsnDraw=False
res2@gsnFrame=False
res2@gsnLeftString=""
res2@gsnRightString=""
res2@tiXAxisString=""
res2@tiYAxisString=""
res2@cnLineLabelsOn=False
res2@cnLevelSpacingF=.05
res2@gsnContourZeroLineThicknessF=0
res2@gsnContourNegLineDashPattern=1
res2@trXMinF=$XMIN
res2@trXMaxF=$XMAX
plot2=gsn_csm_contour(wks,z2,res2)
;plot2=ColorShadeLeGeContour(plot2,-0.01,"cyan",0.01,"salmon")

overlay(plot1,plot2)

draw(plot1)
frame(wks)
end
EOF
ncl tmp.ncl
rm tmp.ncl
mv out.eps figs/bv$BVEXP.eps
display figs/bv$BVEXP.eps
echo "NORMAL END"
