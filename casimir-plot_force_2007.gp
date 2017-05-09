set terminal postscript landscape color enhanced "Times-Roman" 16
set output 'gold_force_decca_2007_dif.eps'
set autoscale 
unset logscale xy
set xtic auto                      
set ytic auto
unset logscale xy
set xrange [150:750]
set yrange [-10:105]
set xlabel "a, nm" font "Times-Roman, 20"
set ylabel "|{P^{th}-P^{exp}|}, mPa" font "Times-Roman, 20" 
set key right top	#set the legend
# unset key

set style line 1 lt 1 lc rgb "#00eeee" lw 2 pt 11 ps 1
set style line 2 lt 1 lc rgb "#008040" lw 2 pt 2 ps 1
set style line 3 lt 1 lc rgb "#ffa040" lw 2 pt 5 ps 1
set style line 4 lt 1 lc rgb "#9400d3" lw 2 pt 2 ps 1
set style line 5 lt 1 lc rgb "#191970" lw 2 pt 7 ps 0.8
set style line 6 lt 0 lc rgb "#ff0000" lw 8 pt 2 ps 1.5
set style line 7 lt 1 lc rgb "#f055f0" lw 2 pt 2 ps 1
set style line 8 lt 0 lc rgb "#804014" lw 7 pt 5 ps 1
set style line 9 lt 1 lc rgb "#000000" lw 3 pt 5 ps 1

#plot "decca-2007-paper-data.txt" u 1:3 w p ls 1 t "Generalized Plasma (paper)",\
#"casimir_decca_2007.txt" u 1:(-$5) smooth csplines \
#w l ls 2 t "Generalized Plasma",\
#"decca-2007-paper-data.txt" u 1:5 w p t "Drude (paper)" ls 5,\
#"casimir_decca_2007.txt" u 1:(-$3) smooth csplines \
#w l ls 3 t "Drude",\
#"decca-2007-paper-data.txt" u 1:2 w p ls 4 t "Experiment (paper)",\
#"casimir_decca_2007.txt" u 1:(-$2) smooth csplines \
#w l ls 6 t "Kramers-Kronig", \
#"casimir_decca_2007.txt" u 1:(-$8) smooth csplines \
#w l ls 7 t "Gauss"

max(x,y) = (x > y) ? x : y

plot "casimir_decca_2007.txt" u 1:(abs($2))  smooth csplines \
w l ls 1 t "Kramers-Kronig",\
"casimir_decca_2007.txt" u 1:(abs($5))  smooth csplines \
w l ls 2 t "Generalized plasma", \
"casimir_decca_2007.txt" u 1:(abs($6)) smooth csplines \
w l ls 3 t "Drude-Lorentz",\
"casimir_decca_2007.txt" u 1:(abs($7))  smooth csplines \
w l ls 4 t "Brendel-Bormann", \
"casimir_decca_2007.txt" u 1:(abs($8))  smooth csplines \
w l  ls 5 t "Gauss", \
"casimir_decca_2007.txt" u 1:(abs($4))  smooth csplines \
w l  ls 6 t "Marachevsky", \
"plot_decca_2007_error_int.csv" u 1:(max(abs($2),abs($3))) smooth csplines \
w l ls 9 t "95% confidence interval"