gnuplot -pe << EOF
set term postscript enhanced
set xlabel "x"
set ylabel "y"
set zlabel "z"
set output 'threedplot.pdf'
set title 'Particle Traveling in 3d' 
splot 'output.dat' using 2:3:4 with lines not
EOF

gnuplot -pe << EOF
set term postscript enhanced
set xlab "time [s]"
set ylab "Velcoty (magnitude) [m/s]"
set output 'velocityvstime.pdf'
set title 'Velocity over time'
plot 'output.dat' using 1:15 with lines not
EOF

gnuplot -pe << EOF
set term postscript enhanced
set xlab "time [s]"
set ylab "Energy loss"
set output 'elossperstep.pdf'
set title 'Energy Loss per Step'
plot 'output.dat' using 1:16 with lines not
EOF

gnuplot -pe << EOF
set xlab "time [s]"
set ylab "Magnetic Field"
set title 'Magnetic field per time'
set st d l
plot "output.dat" u 1:11 tit "Bx", \
     "output.dat" u 1:12 tit "By", \
     "output.dat" u 1:13 tit "Bz", \
     "output.dat" u 1:14 tit "Magntitude of Bfield"
set term pdf enhanced
set output 'Bfieldovertime.pdf'
rep
EOF

#rep 'output.dat' u 1:(sqrt(($11*$11+$12*$12+$13*$13))) tit "Magntitude of Bfield"

#display *.pdf
