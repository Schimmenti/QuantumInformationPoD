set xlabel "Matrix dimension"
set ylabel "Log execution time [s]"
array fnames[7]
fnames[1]="modeijk"
fnames[2]="modeikj"
fnames[3]="modejik"
fnames[4]="modejki"
fnames[5]="modekij"
fnames[6]="modekji"
fnames[7]="intrinsic"
plot for [i=1:7] fnames[i] using 1:2 title fnames[i] with linespoints lw 2
set logscale y
set term png size 1024, 720
set output "time_plots.png"
replot
set term x11
