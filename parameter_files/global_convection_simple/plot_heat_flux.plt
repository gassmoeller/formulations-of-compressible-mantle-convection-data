# Gnuplot script

set title "Heat flux"
set xlabel "Time in years"
set ylabel "Heat flux through top boundary in W/m"
set key bottom samplen 0.0
set format y "%.1tx10^{%T}"
set xrange [1.5e8:4e8]

set terminal postscript color portrait dashed enhanced 'Arial'
set output 'heat_flux.eps'
set size 1.0,0.6

set label "time lag" at 2.27e8,6.5e5 textcolor "#e56b5d"

plot "ALA/output-anelastic-liquid-approximation/statistics" using 2:23 with lines linecolor "#4b03a1" lw 4 title "Anelastic Liquid Approximation", \
	"isentropic/output-isothermal-compression/statistics" using 2:23 with lines dashtype 1 linecolor rgb "#fdc328" lw 4 title "Isentropic compression", \
	"projected-density/output-projected-density-compression/statistics" using 2:23 with lines dashtype 1 linecolor "#e56b5d" lw 4 title "Projected density"
#replot
