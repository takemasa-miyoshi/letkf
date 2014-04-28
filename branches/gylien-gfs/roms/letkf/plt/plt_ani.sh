#!/bin/sh
set -e
END=120
VAR=temp
XMIN=-770
XMAX=0
Z=40
CDIR=`pwd`
cd ../..
ROMS=`pwd`
DATADIR=$ROMS/DATA/nature
OUTNAME=nature_$VAR
TITLE="NATURE SST"
Z1=`expr $Z - 1` || test 1 -eq 1
cd $CDIR
cat << EOF > tmp.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;wks=gsn_open_wks("x11","gsun02n")
wks=gsn_open_wks("ps","out")
;gsn_define_colormap(wks,"testcmap")
gsn_define_colormap(wks,"rainbow")
;gsn_define_colormap(wks,"gui_default")
;gsn_define_colormap(wks,"gsltod")
files=systemfunc("ls $DATADIR/*.nc")
data=addfiles(files,"r")
z=addfiles_GetVar(data,files,"$VAR")
xlab=fspan(-770,21,114)
ylab=fspan(0,791,114)
z&eta_rho=ylab
z&xi_rho=xlab

res1=True
res1@tiMainString="NATURE(shade), BV(contour)"
res1@gsnDraw=False
res1@gsnFrame=False
res1@gsnSpreadColors=True
res1@cnFillOn=True
res1@cnLinesOn=False
res1@cnLevelSelectionMode="ManualLevels"
;res1@cnMaxLevelValF=2.45
;res1@cnMinLevelValF=2.15
;res1@cnLevelSpacingF=0.005
res1@cnMaxLevelValF=20
res1@cnMinLevelValF=10
res1@cnLevelSpacingF=0.1
res1@cnLineLabelsOn=False
res1@lbLabelAutoStride=True
res1@lbOrientation="Vertical"
res1@trXMinF=$XMIN
res1@trXMaxF=$XMAX
plot1=gsn_csm_contour(wks,z(0,$Z1,:,:),res1)

do i=1,$END
  ii=i-1
  it=i*0.25
  setvalues plot1@data
    "sfDataArray" : z(ii,$Z1,:,:)
  end setvalues
  setvalues plot1
    "tiMainString" : "$TITLE TIME:"+sprintf("%5.2f",it)
  end setvalues
  draw(plot1)
  frame(wks)
end do
end
EOF
ncl tmp.ncl
rm tmp.ncl
convert -delay 10 out.ps out.mpg
#rm out.ps
mv out.mpg figs/$OUTNAME.mpg
echo "NORMAL END"
