#!/bin/sh
#  plot fcst norm
set -e
REC=2
COLOR=color
CDIR=`pwd`
#
cat << EOF > tmp.plt
set key bottom right
set size 0.6
set xtics 360
set xlabel "day"
set ylabel "norm [m/s]"
set title "sqrt(ubar^2 + vbar^2)"
plot "spup_uv.dat" usi $REC notitle w l
set terminal postscript eps enhanced $COLOR
set output "spup_uv.eps"
plot "spup_uv.dat" usi $REC notitle w l
set terminal png
set output "spup_uv.png"
plot "spup_uv.dat" usi $REC notitle w l
EOF
gnuplot -persist tmp.plt
rm tmp.plt
echo "NORMAL END"
