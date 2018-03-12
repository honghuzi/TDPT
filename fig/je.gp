reset

set loadpath 'D:/libs/gnuplotconf'
load 'configs/xyborder.cfg'
# load 'configs/grid.cfg'
load 'palettes/moreland.pal'
#set palette negative
#moreland, parula, orrd, plasma, rdbu, rdylbu, spectral, whylrd, ylrd, bugn, bupu, chromajs, inferno, pubu, purd, rdbu,
set t  eps font "Arial, 14" size 3, 2.8
set size ratio 1
set o "temp"

set title 'joint energy spectrum'
set xlabel "E_1 (eV)"
set ylabel "E_2 (eV)"
unset key
set logscale zcb
# set xrange [0:60]
# set yrange [0:60]

set format cb"10^{%1T}"

set border lw 2
set tics front
# set xtics 1; set ytics 1
# set arrow 1 from -2.5, 0 to 2.5, 0 nohead front ls 101 lw 3
# set arrow 2 from 0, -2.5 to 0, 2.5 nohead front ls 101 lw 3
#set palette defined (0 "blue", 1 "white", 2 "red")
#set view map
plot '../data/je.bin' binary format='%float64' w image
unset output
set o "je.eps"
offset = GPVAL_CB_MAX
# #show variables all
# #set pm3d map
# set cbrange [1e-7: 1e-3]
plot "" binary format='%float64' u 1:2:($3/offset) w image
#splot "" binary format='%float64' u 1:2:($3/offset) w points palette
set output