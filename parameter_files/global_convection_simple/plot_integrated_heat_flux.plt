# Gnuplot script

set title "Heat flux difference integrated over time"
set xlabel "Time in years"
set ylabel "Relative difference in cumulative heat flux"
set key samplen 0.0
set yrange [-0.01:0.02]

set terminal postscript color portrait dashed enhanced 'Arial'
set output 'cumulative_heat_flux.eps'
set size 1.0,0.6

plot "ALA/output-anelastic-liquid-approximation/statistics" using 2:(($23*$3-6.87e5*$3)/6.87e15) with lines smooth cumulative linecolor "#4b03a1" lw 4 title "Anelastic Liquid Approximation", \
	"isentropic/output-isentropic-compression/statistics" using 2:(($23*$3-6.87e5*$3)/6.87e15) with lines smooth cumulative dashtype 1 linecolor rgb "#fdc328" lw 4 title "Isentropic compression", \
	"projected-density/output-projected-density-compression/statistics" using 2:(($23*$3-6.87e5*$3)/6.87e15) with lines smooth cumulative dashtype 1 linecolor "#e56b5d" lw 4 title "Projected density"
#replot
