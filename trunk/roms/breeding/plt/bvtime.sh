#!/bin/sh
#  plot bred vectors
set -e
END=240
MEM=001
XMIN=-770
XMAX=0
Z=40
Y=57
CDIR=`pwd`
cd ../..
ROMS=`pwd`
NATURE=$ROMS/DATA/nature
BVEXP=bv_ubar_5E-3
END1=`expr $END - 1`
Z1=`expr $Z - 1`
Y1=`expr $Y - 1`
cd $CDIR
cat << EOF > tmp.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;wks=gsn_open_wks("x11","gsun02n")
wks=gsn_open_wks("eps","out")
;gsn_define_colormap(wks,"testcmap")
;gsn_define_colormap(wks,"gui_default")
gsn_define_colormap(wks,"rainbow")
;gsn_define_colormap(wks,"gsltod")
files1=systemfunc("ls $NATURE/*.nc")
files2=systemfunc("ls $ROMS/DATA/$BVEXP/bv/$MEM/*.nc")
data1=addfiles(files1,"r")
data2=addfiles(files2,"r")
xlab=fspan(-770,21,114)
ylab=fspan(0,0.25*$END1,$END)

z=addfiles_GetVar(data1,files1,"temp")
z1=z(0:$END1,$Z1,$Y1,:)
z1&xi_rho=xlab
z1&time=ylab
;printVarSummary(z1)
res1=True
res1@tiMainString="NATURE(shade), BV(contour)"
res1@tiYAxisString="TIME [day]"
res1@tiXAxisString="DISTANCE FROM COAST [km]"
res1@gsnDraw=False
res1@gsnFrame=False
res1@gsnSpreadColors=True
res1@cnFillOn=True
res1@cnLinesOn=False
res1@cnLevelSelectionMode="ManualLevels"
res1@cnMaxLevelValF=20
res1@cnMinLevelValF=10
res1@cnLevelSpacingF=0.1
res1@cnLineLabelsOn=False
res1@lbLabelAutoStride=True
res1@lbOrientation="Vertical"
res1@trXMinF=$XMIN
res1@trXMaxF=$XMAX
plot1=gsn_csm_contour(wks,z1,res1)

zz=addfiles_GetVar(data2,files2,"temp")
z2=zz(0:$END1,$Z1,$Y1,:)
z2=z2/max(abs(z2))
z2&xi_rho=xlab
z2&time=ylab
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
res2@cnLevelSpacingF=0.05
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
mv out.eps figs/bvtime$BVEXP.eps
display figs/bvtime$BVEXP.eps
echo "NORMAL END"
