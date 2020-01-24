#!/usr/bin/gnuplot -persist

set terminal postscript eps

reset
set output "xy_trajectory.eps"
set xlabel "x"
set ylabel "y"
set grid
plot 'trace.dat' u 2:3 w l lt 1 lw 2 lc rgb "blue" title "xy trajectory"

reset
set output "xz_trajectory.eps"
set xlabel "x"
set ylabel "z"
set grid
plot 'trace.dat' u 2:4 w l lt 1 lw 2 lc rgb "blue" title "xz trajectory"

reset
set output "velocity.eps"
set xlabel "t"
set ylabel "speed"
set yrange [-1 : 1]
set grid
plot 'trace.dat' u 1:5 w l lt 1 lw 2 lc rgb "blue" title "vx",\
'trace.dat' u 1:6 w l lt 1 lw 2 lc rgb "black" title "vy",\
'trace.dat' u 1:7 w l lt 1 lw 2 lc rgb "red" title "vz"

reset
set output "mass.eps"
set xlabel "t"
set ylabel "mass"
set yrange [0 : 1]
set grid
plot 'trace.dat' u 1:8 w l lt 1 lw 2 lc rgb "blue" title "mass"

reset
set output "control.eps"
set xlabel "t"
set ylabel "u"
set yrange [-1.5 : 1.5]
set grid
plot 'trace.dat' u 1:16 w l lt 1 lw 2 lc rgb "blue" title "ux",\
'trace.dat' u 1:17 w l lt 1 lw 2 lc rgb "black" title "uy",\
'trace.dat' u 1:18 w l lt 1 lw 2 lc rgb "red" title "uz",\
-1 lt 2 lw 2 lc rgb "blue" notitle,\
1 lt 2 lw 2 lc rgb "blue" notitle

reset
set output "velocityNorm.eps"
set xlabel "t"
set ylabel "speed norm"
set grid
plot 'trace.dat' u 1:(sqrt($5*$5+$6*$6+$7*$7)) w l lt 1 lw 2 lc rgb "blue" title "Velocity norm"

reset
set output "hamiltonian.eps"
set xlabel "t"
set ylabel "hamiltonian"
set grid
plot 'trace.dat' u 1:19 w l lt 1 lw 2 lc rgb "blue" title "Hamiltonian"

reset
set output "controlNorm.eps"
set xlabel "t"
set ylabel "u"
set yrange [0 : 1.5]
set grid
plot 'trace.dat' u 1:(sqrt($16*$16+$17*$17+$18*$18)) w l lt 1 lw 2 lc rgb "blue" title "Norm of the control" ,\
1 lt 2 lw 2 lc rgb "blue" notitle

