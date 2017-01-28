 set term png
 set output "./plots/pressure0120.png"
 set xlabel "x"
 set ylabel "pressure"
 m="./data/frame0120t_0.0999.dat"
 set nokey
 set grid
 set title 'Frame 0120, time = 0.0999'
 plot m using 1:0003 with linespoints
