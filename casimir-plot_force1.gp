set terminal postscript landscape color enhanced "Times-Roman" 16
set output 'gold_force_decca.eps'
set autoscale 
unset logscale xy
set xtic auto                      
set ytic auto
set xrange [200:700]
set xlabel "a, nm" font "Times-Roman, 20"
set ylabel "Force, fN" font "Times-Roman, 20" 
set key right bottom	#set the legend
# unset key

set style line 1 lt 1 lc rgb "#00eeee" lw 2.5 pt 7 ps 1
set style line 2 lt 1 lc rgb "#008040" lw 1.5 pt 2 ps 1
set style line 3 lt 1 lc rgb "#ffa040" lw 2 pt 5 ps 1
set style line 4 lt 1 lc rgb "#9400d3" lw 2 pt 5 ps 1
set style line 5 lt 0 lc rgb "#191970" lw 8 pt 7 ps 1
set style line 6 lt 0 lc rgb "#ff0000" lw 8 pt 2 ps 1.5
set style line 7 lt 0 lc rgb "#f055f0" lw 7 pt 2 ps 1
set style line 8 lt 0 lc rgb "#804014" lw 7 pt 5 ps 1

plot "casimir_decca.txt" u 1:2 smooth csplines w l ls 1 t "Drude",\
"casimir_decca.txt" u 1:3 smooth csplines \
w l ls 2 t "Plasma",\
"Decca-2016-06-au-37.csv" u 1:2 w p t "Experiment"