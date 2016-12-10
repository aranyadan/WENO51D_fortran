set xlabel "x"
set ylabel "y"
m="./data"
set terminal x11 0
set nokey
set grid
set title 'IC'
plot m using 1:2 with linespoints
set term png
set output "./plots/out.png"
replot
