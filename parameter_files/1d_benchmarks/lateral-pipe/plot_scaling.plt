# Gnuplot script

set terminal pdf color solid dashed font "Arial,12" size 30cm,10cm
set output 'mass_flux_error.pdf'

set multiplot layout 1,3
set size ratio 1.0


set title "Mass flux error - Adiabatic temperature"
set xlabel "h"
set ylabel "Mass flux error"
set logscale xy
set format y "10^{%S}"
set ytics (1,1e-3,1e-6,1e-9)

set xrange [0.5:100]
set yrange [1e-9:10]

set grid ytics
set key top right noautotitles

f(x)= 2e-4 / (x*x)

plot "mass_flux_error_adiabatic" using 1:2 with linespoints linetype 1 pointtype 7 ps 0.7 linecolor rgb "#0A5F02" lw 3 title "Anelastic liquid approximation (ALA)", \
     "mass_flux_error_adiabatic" using 1:3 with linespoints linetype 1 pointtype 9 ps 0.7 linecolor rgb "#1c567a" lw 3 title "Isothermal compression", \
     "mass_flux_error_adiabatic" using 1:4 with linespoints linetype 1 pointtype 11 ps 0.7 linecolor rgb "#814292" lw 3 title "Hydrostatic compression", \
     "mass_flux_error_adiabatic" using 1:5 with linespoints linetype 1 pointtype 13 ps 0.7 linecolor rgb "#d7383b" lw 3 title "Projected density", \
     f(x) title '1/h^2' with lines linestyle 2 linecolor "gray" lw 3 



set size ratio 1.0
#set origin 0.0,0.5
set key off

set title "Mass flux error - Sub-adiabatic temperature"
set xlabel "h"
set ylabel "Mass flux error"
set logscale xy

set xrange [0.5:100]
set yrange [1e-9:10]

set grid ytics

set key top right noautotitles

plot "mass_flux_error_sub-adiabatic" using 1:2 with linespoints linetype 4 pointtype 7 ps 0.7 linecolor rgb "#0A5F02" lw 3 title "Anelastic liquid approximation (ALA)", \
     "mass_flux_error_sub-adiabatic" using 1:3 with linespoints linetype 4 pointtype 9 ps 0.7 linecolor rgb "#1c567a" lw 3 title "Isothermal compression", \
     "mass_flux_error_sub-adiabatic" using 1:4 with linespoints linetype 4 pointtype 11 ps 0.7 linecolor rgb "#814292" lw 3 title "Hydrostatic compression", \
     "mass_flux_error_sub-adiabatic" using 1:5 with linespoints linetype 4 pointtype 13 ps 0.7 linecolor rgb "#d7383b" lw 3 title "Projected density", \
     f(x) title '1/h^2' with lines linestyle 2 linecolor "gray" lw 3 

set size ratio 1.0
#set origin 0.0,0.5
set key off

set title "Mass flux error - Super-adiabatic temperature"
set xlabel "h"
set ylabel "Mass flux error"
set logscale xy

set xrange [0.5:100]
set yrange [1e-9:10]

set grid ytics

set key top right noautotitles

plot "mass_flux_error_super-adiabatic" using 1:2 with linespoints linetype 4 pointtype 7 ps 0.7 linecolor rgb "#0A5F02" lw 3 title "Anelastic liquid approximation (ALA)", \
     "mass_flux_error_super-adiabatic" using 1:3 with linespoints linetype 4 pointtype 9 ps 0.7 linecolor rgb "#1c567a" lw 3 title "Isothermal compression", \
     "mass_flux_error_super-adiabatic" using 1:4 with linespoints linetype 4 pointtype 11 ps 0.7 linecolor rgb "#814292" lw 3 title "Hydrostatic compression", \
     "mass_flux_error_super-adiabatic" using 1:5 with linespoints linetype 4 pointtype 13 ps 0.7 linecolor rgb "#d7383b" lw 3 title "Projected density", \
     f(x) title '1/h^2' with lines linestyle 2 linecolor "gray" lw 3 

unset multiplot
#replot
