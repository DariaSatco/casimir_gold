set terminal postscript landscape color enhanced "Times-Roman" 16
set output 'gold_force_lamoreaux.eps'
set autoscale 
set logscale x
set xtic auto                      
set ytic auto
set xrange [0.7:7]
set xtics(0.7,0.8,0.9,1,2,3,4,5,6,7)
set xlabel "a, {/Symbol m}m" font "Times-Roman, 20"
set ylabel "Force, pN*{/Symbol m}m^2" font "Times-Roman, 20" 
set key right bottom	#set the legend
# unset key

set style line 1 lt 1 lc rgb "#00eeee" lw 2.5 pt 7 ps 1
set style line 2 lt 1 lc rgb "#008040" lw 1.5 pt 2 ps 1
set style line 3 lt 1 lc rgb "#ffa040" lw 2 pt 5 ps 1
set style line 4 lt 1 lc rgb "#9400d3" lw 2 pt 5 ps 1
set style line 5 lt 0 lc rgb "#191970" lw 8 pt 7 ps 1
set style line 6 lt 0 lc rgb "#ff0000" lw 8 pt 7 ps 1.5
set style line 7 lt 1 lc rgb "#f055f0" lw 2 pt 2 ps 1
set style line 8 lt 0 lc rgb "#804014" lw 7 pt 5 ps 1

plot "casimir_lamoreaux.txt" u 1:2 smooth csplines w l ls 1 t "Kramers-Kronig",\
"casimir_lamoreaux.txt" u 1:3 smooth csplines \
w l ls 2 t "Drude",\
"casimir_lamoreaux.txt" u 1:4 smooth csplines \
w l ls 3 t "Marachevsky",\
"casimir_lamoreaux.txt" u 1:5 smooth csplines \
w l ls 4 t "Generalized Plasma",\
"casimir_lamoreaux.txt" u 1:6 smooth csplines \
w l ls 5 t "Drude-Lorentz",\
"casimir_lamoreaux.txt" u 1:7 smooth csplines \
w l ls 6 t "Brendel-Bormann",\
"casimir_lamoreaux.txt" u 1:8 smooth csplines \
w l ls 7 t "Gauss",\
"lamoreaux-2010-fig2-1.csv" u 1:2:3:4 with errorbars t "Experiment"